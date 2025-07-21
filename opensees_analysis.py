import numpy as np
import zipfile
import xml.etree.ElementTree as ET
try:
    from openseespy import opensees as ops
except Exception:
    ops = None

def _require_ops():
    if ops is None:
        raise RuntimeError(
            'OpenSeesPy is required for structural analysis. '
            'Install it with "pip install openseespy" and rerun.')

# section tables duplicated from anly.py
list_H = [400,450,500,550,600,650,700,750,800,850,900,950,1000]
list_B = [200,250,300,350,400]
list_tw = [9,12,14,16,19,22]
list_tf = [12,16,19,22,25,28,32,36,40]
list_D = [350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
list_t = [12,16,19,22,25,28,32,36,38,40]
E = 2.05e5
G = E/(2*(1+0.3))
FY_BEAM = 325.0
FY_COL = 235.0

def _read_params(path='param_3d.xlsx'):
    vals = {}
    with zipfile.ZipFile(path) as zf:
        data = zf.read('xl/worksheets/sheet1.xml')
    root = ET.fromstring(data)
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    for c in root.findall('.//m:c', ns):
        r = c.attrib.get('r')
        v = c.find('m:v', ns)
        if v is not None:
            vals[r] = float(v.text)
    def cell(addr, default=0.0):
        return vals.get(addr, default)

    xSpan = int(cell('A1', 1))
    ySpan = int(cell('B1', xSpan))
    zSpan = int(cell('C1', 1))
    xlen  = cell('D1', 6000.0)
    ylen  = cell('E1', 6000.0)
    story_def = cell('F1', 4000.0)
    load  = cell('G1', 0.01)  # already in N/mm^2
    param = cell('H1', 0.0)
    heights = []
    # row 2 contains individual story heights A2,B2,...
    def col_letter(idx):
        s=""
        while True:
            s = chr(ord('A') + idx % 26) + s
            idx //= 26
            if idx==0:
                break
            idx -= 1
        return s
    for i in range(zSpan):
        caddr = f"{col_letter(i)}2"
        heights.append(cell(caddr, story_def))
    return xSpan, ySpan, zSpan, xlen, ylen, heights, load, param

def _code2idz(code, cb):
    if code == '0':
        return [0,0,0,0]
    if cb == 0:
        pos = str(code).zfill(2)
        return [int(pos[0],16), int(pos[1],16),0,0]
    else:
        pos = str(code).zfill(4)
        return [int(pos[0],16), int(pos[1],16), int(pos[2],16), int(pos[3],16)]

def _col_props(code, col_list):
    cid = str(code)[2:]
    if cid in col_list:
        d = col_list[cid]
        D = float(d['D'])
        td = float(d['td'])
        A = float(d['A'])
        Z = float(d['Z'])
    else:
        vals = _code2idz(cid,0)
        D = list_D[vals[0]-1]
        td = list_t[vals[1]-1]
        A = D*D - (D-2*td)**2
        Z = (D**4 - (D-2*td)**4)/(6*D)
    Iy = Z*(D/2)
    return {'A':A,'Iy':Iy,'Iz':Iy,'J':2*Iy,'Z':Z,'fy':FY_COL}

def _beam_props(code, beam_list):
    bid = str(code)[2:]
    if bid in beam_list:
        d = beam_list[bid]
        H = float(d['H'])
        B = float(d['B'])
        tw = float(d['tw'])
        tf = float(d['tf'])
        A = float(d['A'])
        Z = float(d['Z'])
    else:
        vals = _code2idz(bid,1)
        H = list_H[vals[0]-1]
        B = list_B[vals[1]-1]
        tw = list_tw[vals[2]-1]
        tf = list_tf[vals[3]-1]
        A = H*B - (B-tw)*(H-2*tf)
        Z = (B*H**3 - (B-tw)*(H-2*tf)**3)/(6*H)
    Iy = Z*(H/2)
    Iz = Z*(B/2)
    return {'A':A,'Iy':Iy,'Iz':Iz,'J':Iy+Iz,'Z':Z,'fy':FY_BEAM}

def _stress_ratio(axial, moment, prop):
    sigma_a = abs(axial)/prop['A']
    sigma_b = abs(moment)/prop['Z']
    return max(sigma_a, sigma_b)/prop['fy']

def run_analysis(member, col_list, beam_list, arraymat, xSpan, ySpan, zSpan, margin=2):
    _require_ops()
    p = _read_params()
    xSpan, ySpan, zSpan, xlen, ylen, heights, load, _ = p
    lx = 2*int(xSpan)+1
    ly = 2*int(ySpan)+1
    ops.wipe()
    ops.model('basic','-ndm',3,'-ndf',6)
    # use separate orientation vectors to avoid parallel configurations
    # transfTag 1: columns and beams along the global X or Z axis
    #   orientation vector is (0, 1, 0) so it is never parallel to the element
    ops.geomTransf('Linear', 1, 0, 1, 0)
    # transfTag 2 is for beams along the global Y-axis
    ops.geomTransf('Linear', 2, 1, 0, 0)
    node_map = {}
    tag=1
    cumh = [0.0]
    for h in heights:
        cumh.append(cumh[-1] + h)
    for iz in range(zSpan+1):
        zc = cumh[iz]
        for jy in range(0,ly,2):
            for kx in range(0,lx,2):
                x = (kx//2)*xlen
                y = (jy//2)*ylen
                ops.node(tag,x,y,zc)
                if iz==0:
                    ops.fix(tag,1,1,1,1,1,1)
                node_map[(iz,jy,kx)] = tag
                tag+=1
    eleTag=1
    bx_tags=[]; by_tags=[]; c_tags=[]
    bx_props=[]; by_props=[]; c_props=[]
    bx_story={}
    by_story={}
    # columns
    for iz in range(zSpan):
        for jy in range(0,ly,2):
            for kx in range(0,lx,2):
                code = member[iz+margin, jy+margin, kx+margin]
                prop = _col_props(code, col_list)
                n1=node_map[(iz,jy,kx)]
                n2=node_map[(iz+1,jy,kx)]
                ops.element('elasticBeamColumn',eleTag,n1,n2,prop['A'],E,G,prop['J'],prop['Iy'],prop['Iz'],1)
                c_tags.append(eleTag)
                c_props.append(prop)
                eleTag+=1
    # beams x
    for iz in range(zSpan):
        for jy in range(0,ly,2):
            for kx in range(1,lx-1,2):
                code = member[iz+margin, jy+margin, kx+margin]
                prop = _beam_props(code, beam_list)
                n1=node_map[(iz,jy,kx-1)]
                n2=node_map[(iz,jy,kx+1)]
                ops.element('elasticBeamColumn',eleTag,n1,n2,prop['A'],E,G,prop['J'],prop['Iy'],prop['Iz'],1)
                bx_tags.append(eleTag)
                bx_props.append(prop)
                bx_story.setdefault(iz+1, []).append(eleTag)
                eleTag+=1
    # beams y
    for iz in range(zSpan):
        for jy in range(1,ly-1,2):
            for kx in range(0,lx,2):
                code = member[iz+margin, jy+margin, kx+margin]
                prop = _beam_props(code, beam_list)
                n1=node_map[(iz,jy-1,kx)]
                n2=node_map[(iz,jy+1,kx)]
                ops.element('elasticBeamColumn',eleTag,n1,n2,prop['A'],E,G,prop['J'],prop['Iy'],prop['Iz'],2)
                by_tags.append(eleTag)
                by_props.append(prop)
                by_story.setdefault(iz+1, []).append(eleTag)
                eleTag+=1
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    w_x = load * ylen
    w_y = load * xlen
    for iz,tags in bx_story.items():
        for tag in tags:
            ops.eleLoad('-ele', tag, '-type', '-beamUniform', 0, 0, -w_x)
    for iz,tags in by_story.items():
        for tag in tags:
            ops.eleLoad('-ele', tag, '-type', '-beamUniform', 0, 0, -w_y)
    # show the uniform line load values applied along each beam direction
    print(f"Applied load {load} N/mm^2 -> w_x={w_x:.2f} N/mm, w_y={w_y:.2f} N/mm")
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('LoadControl',1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)
    deflect=[]
    for iz in range(1,zSpan+1):
        vals=[]
        for jy in range(0,ly,2):
            for kx in range(0,lx,2):
                vals.append(abs(ops.nodeDisp(node_map[(iz,jy,kx)],3)))
        deflect.append(max(vals) if vals else 0.0)
    stress_bx=[]
    for tag,prop in zip(bx_tags,bx_props):
        f=ops.eleForce(tag)
        m=max(abs(f[4]),abs(f[5]),abs(f[10]),abs(f[11]))
        stress_bx.append(_stress_ratio(f[0],m,prop))
    stress_by=[]
    for tag,prop in zip(by_tags,by_props):
        f=ops.eleForce(tag)
        m=max(abs(f[4]),abs(f[5]),abs(f[10]),abs(f[11]))
        stress_by.append(_stress_ratio(f[0],m,prop))
    stress_c=[]; cof=[]
    for tag,prop in zip(c_tags,c_props):
        f=ops.eleForce(tag)
        m=max(abs(f[4]),abs(f[5]),abs(f[10]),abs(f[11]))
        sr=_stress_ratio(f[0],m,prop)
        stress_c.append(sr)
        ratio = abs(f[0])/(prop['A']*prop['fy'])
        cof.append([ratio, ratio])
    cof=np.array(cof).T
    return np.array(stress_bx), np.array(stress_by), np.array(stress_c), cof, np.array(deflect)