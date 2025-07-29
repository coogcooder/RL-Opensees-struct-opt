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


def _read_params(path='csv/param_3d.xlsx'):
    """Read parameters from Excel file - exactly like MATLAB version"""
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
    # plastic section modulus
    zpy = D*td*(D-td) + (D-2*td)**2*td/2
    return {'A':A,'Iy':Iy,'Iz':Iy,'J':2*Iy,'Z':Z,'fy':FY_COL,'zpy':zpy}


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
    zpy = B*tf*(H-tf) + (H-2*tf)**2*tw/4
    return {'A':A,'Iy':Iy,'Iz':Iz,'J':Iy+Iz,'Z':Z,'fy':FY_BEAM,'zpy':zpy}


def _stress_ratio(axial, moment, prop):
    sigma_a = abs(axial)/prop['A']
    sigma_b = abs(moment)/prop['Z']
    return max(sigma_a, sigma_b)/prop['fy']


def _earthquake(zSpan, pw, heights, x_total, y_total, Z=1.0, Co=0.2):
    """Return lateral seismic forces for each story."""
    T = 0.03 * sum(heights) * 1e-3
    wi = x_total * y_total * (pw - 1e-3)
    if T < 0.6:
        Rt = 1.0
    elif 0.6 < T < 2 * 0.6:
        Rt = 1 - 0.2 * (T / 0.6 - 1) ** 2
    else:
        Rt = 1.6 * 0.6 / T

    Sigmawi = [0.0] * zSpan
    for i in range(zSpan - 1, -1, -1):
        if i == zSpan - 1:
            Sigmawi[i] = wi
        else:
            Sigmawi[i] = Sigmawi[i + 1] + wi

    Qi = [0.0] * zSpan
    for i in range(zSpan):
        alpha = Sigmawi[i] / Sigmawi[0]
        Ai = 1 + (1 / np.sqrt(alpha) - alpha) * (2 * T / (1 + 3 * T))
        Ci = Z * Rt * Ai * Co
        Qi[i] = Ci * Sigmawi[i]

    seismic = [0.0] * zSpan
    for i in range(zSpan - 1, -1, -1):
        if i == zSpan - 1:
            seismic[i] = Qi[i]
        else:
            seismic[i] = Qi[i] - Qi[i + 1]
    return np.array(seismic)


def _proof_stress(zpy, js, je, jel, c_g, frame, nelx, nely, nelz, nj, ng, F=None):
    """Replicate MATLAB ``proof_stress`` behavior."""
    if F is None:
        F = [FY_COL, FY_BEAM]
    zpy = np.asarray(zpy, dtype=float)
    js = np.asarray(js, dtype=int)
    je = np.asarray(je, dtype=int)
    jel = np.asarray(jel, dtype=int)
    c_g = np.asarray(c_g, dtype=int)

    if frame.startswith('s'):
        nfj = (nelx + 1) * (nely + 1)
        ngx = nelx * (nely + 1) * nelz
    elif frame.startswith('x'):
        nfj = nelx + 1
        ngx = ng
    else:
        nfj = nely + 1
        ngx = 0

    nrps = nj - nfj * 2
    rps = np.zeros(nrps * 2)

    for i in range(nfj + 1, nfj + nrps + 1):
        men = np.where((js == i) | (je == i))[0]
        pci = []
        pbix = []
        pbiy = []
        for m in men:
            strength = F[jel[m] - 1]
            if c_g[m] == 1:
                pci.append(zpy[m] * strength)
            elif c_g[m] == 2 and m < ngx:
                pbix.append(zpy[m] * strength)
            else:
                pbiy.append(zpy[m] * strength)
        rpc = float(np.sum(pci))
        rpbx = float(np.sum(pbix))
        rpby = float(np.sum(pbiy))
        idx = (i - nfj - 1) * 2
        if rpc != 0:
            rps[idx] = 1.5 * rpbx / rpc - 1
            rps[idx + 1] = 1.5 * rpby / rpc - 1
    return rps


def variable_to_member(variable_, xSpan, ySpan, zSpan, margin=2):
    """Convert MATLAB-style variable array to member array format"""
    variable_ = np.array(variable_).flatten()
    
    # Calculate array dimensions
    lx = 2 * xSpan + 1
    ly = 2 * ySpan + 1
    
    # Create arraymat (same as MATLAB framemat)
    arraymat = np.ones((zSpan, ly, lx)) * 2
    for i in range(zSpan):
        for j in range(ly):
            for k in range(lx):
                if j % 2 == 0 and k % 2 == 0:
                    arraymat[i, j, k] = 1
                elif j % 2 == 1 and k % 2 == 1:
                    arraymat[i, j, k] = 0
    
    # Create member array with margins
    member = np.full((zSpan + 2*margin, ly + 2*margin, lx + 2*margin), '0', dtype=object)
    
    # Calculate how many beams and columns we have
    ng = 0  # number of beams
    nc = 0  # number of columns
    for i in range(zSpan):
        for j in range(ly):
            for k in range(lx):
                if arraymat[i, j, k] == 2:  # beam
                    ng += 1
                elif arraymat[i, j, k] == 1:  # column
                    nc += 1
    
    # Extract variables (H, B, tw, tf for beams; D, t for columns)
    nvg = ng  # number of beam groups
    nvc = nc  # number of column groups
    
    # Beam variables: H, B, tw, tf
    H_vals = variable_[:nvg]
    B_vals = variable_[nvg:2*nvg]
    tw_vals = variable_[2*nvg:3*nvg]
    tf_vals = variable_[3*nvg:4*nvg]
    
    # Column variables: D, t
    D_vals = variable_[4*nvg:4*nvg+nvc]
    t_vals = variable_[4*nvg+nvc:4*nvg+2*nvc]
    
    # Convert to discrete indices
    def continuous_to_discrete(val, val_list):
        idx = int(np.clip(np.round(val), 1, len(val_list)))
        return idx
    
    # Fill member array
    beam_idx = 0
    col_idx = 0
    
    for i in range(zSpan):
        for j in range(ly):
            for k in range(lx):
                if arraymat[i, j, k] == 2:  # beam
                    H_idx = continuous_to_discrete(H_vals[beam_idx], list_H)
                    B_idx = continuous_to_discrete(B_vals[beam_idx], list_B)
                    tw_idx = continuous_to_discrete(tw_vals[beam_idx], list_tw)
                    tf_idx = continuous_to_discrete(tf_vals[beam_idx], list_tf)
                    
                    code = hex(H_idx)[2:] + hex(B_idx)[2:] + hex(tw_idx)[2:] + hex(tf_idx)[2:]
                    member[i + margin, j + margin, k + margin] = '0x' + code
                    beam_idx += 1
                    
                elif arraymat[i, j, k] == 1:  # column
                    D_idx = continuous_to_discrete(D_vals[col_idx], list_D)
                    t_idx = continuous_to_discrete(t_vals[col_idx], list_t)
                    
                    code = hex(D_idx)[2:] + hex(t_idx)[2:]
                    member[i + margin, j + margin, k + margin] = '0x' + code
                    col_idx += 1
    
    return member, arraymat


def run_analysis(member, col_list, beam_list, arraymat, xSpan, ySpan, zSpan, margin=2):
    """Main analysis function - EXACTLY matches MATLAB execute_stkcv output format"""
    _require_ops()
    
    # Read parameters
    try:
        p = _read_params()
        xSpan, ySpan, zSpan, xlen, ylen, heights, load, _ = p
    except:
        # Fallback values if file reading fails
        xlen, ylen = 6000.0, 6000.0
        heights = [4000.0] * zSpan
        load = 0.01
    
    lx = 2*int(xSpan)+1
    ly = 2*int(ySpan)+1
    
    # Initialize OpenSees
    ops.wipe()
    ops.model('basic','-ndm',3,'-ndf',6)
    
    # Define transformations
    ops.geomTransf('Linear', 1, 0, 1, 0)
    ops.geomTransf('Linear', 2, 1, 0, 0)
    
    # Create nodes
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
    
    # Create elements
    eleTag=1
    bx_tags=[]; by_tags=[]; c_tags=[]
    bx_props=[]; by_props=[]; c_props=[]
    bx_story={}; by_story={}
    js=[]; je=[]; jel=[]; c_g=[]; zpy=[]
    
    # Columns
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
                js.append(n1)
                je.append(n2)
                jel.append(1)
                c_g.append(1)
                zpy.append(prop['zpy'])
                eleTag+=1
    
    # Beams X-direction
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
                js.append(n1)
                je.append(n2)
                jel.append(2)
                c_g.append(2)
                zpy.append(prop['zpy'])
                eleTag+=1
    
    # Beams Y-direction
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
                js.append(n1)
                je.append(n2)
                jel.append(2)
                c_g.append(2)
                zpy.append(prop['zpy'])
                eleTag+=1
    
    # Apply loads
    ops.timeSeries('Linear',1)
    ops.pattern('Plain',1,1)
    
    # Distributed loads on beams
    w_x = load * ylen
    w_y = load * xlen
    for tags in bx_story.values():
        for tag in tags:
            ops.eleLoad('-ele', tag, '-type', '-beamUniform', 0, 0, -w_x)
    for tags in by_story.values():
        for tag in tags:
            ops.eleLoad('-ele', tag, '-type', '-beamUniform', 0, 0, -w_y)

    # Seismic loads
    seismic = _earthquake(zSpan, load, heights, xSpan * xlen, ySpan * ylen)
    nodes_per_floor = (xSpan + 1) * (ySpan + 1)
    for iz in range(1, zSpan + 1):
        fx = seismic[iz - 1] / nodes_per_floor
        for jy in range(0, ly, 2):
            for kx in range(0, lx, 2):
                ops.load(node_map[(iz, jy, kx)], fx, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Run analysis
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('LoadControl',1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ops.analyze(1)
    
    # Extract results
    # Deflections
    deflect=[]
    for iz in range(1,zSpan+1):
        vals=[]
        for jy in range(0,ly,2):
            for kx in range(0,lx,2):
                vals.append(abs(ops.nodeDisp(node_map[(iz,jy,kx)],3)))
        deflect.append(max(vals) if vals else 0.0)
    
    # Beam stresses
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
    
    # Column stresses and COF
    stress_c=[]; cof=[]
    for tag,prop in zip(c_tags,c_props):
        f=ops.eleForce(tag)
        m=max(abs(f[4]),abs(f[5]),abs(f[10]),abs(f[11]))
        sr=_stress_ratio(f[0],m,prop)
        stress_c.append(sr)
        ratio = abs(f[0])/(prop['A']*prop['fy'])
        cof.append([ratio, ratio])
    
    cof=np.array(cof).T
    
    # Clean up
    ops.wipe()
    
    # Return in EXACT MATLAB format
    return (
        np.array(stress_bx),
        np.array(stress_by), 
        np.array(stress_c),
        cof,
        np.array(deflect).reshape(1, -1)  # Match MATLAB 2D format
    )


def execute_stkcv(choice_section):
    """Direct replacement for MATLAB execute_stkcv function"""
    from module import envmnt_lrn
    
    # Get current environment parameters
    xSpan, ySpan, zSpan, xlen, ylen, heights, load, _ = _read_params()
    
    # Convert choice_section to member format
    member, arraymat = variable_to_member(choice_section, xSpan, ySpan, zSpan)
    
    # Get material lists
    col_list = envmnt_lrn.col_list
    beam_list = envmnt_lrn.beam_list
    
    # Run analysis
    stress_bx, stress_by, stress_c, cof, deflect = run_analysis(
        member, col_list, beam_list, arraymat, xSpan, ySpan, zSpan
    )
    
    return stress_bx, stress_by, stress_c, cof, deflect
