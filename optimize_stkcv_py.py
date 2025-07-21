import numpy as np
import opensees_analysis as osa

def optimize_stkcv(member, col_list, beam_list, arraymat, xSpan, ySpan, zSpan, margin=2):
    """Python replacement for MATLAB optimize_stkcv.

    For now this wrapper simply runs a single OpenSees analysis using
    ``osa.run_analysis`` and returns the resulting stress ratios, capacity
    factors and deflections. It provides a consistent interface so the
    reinforcement learning environment remains unchanged.
    """
    return osa.run_analysis(member, col_list, beam_list, arraymat,
                             xSpan, ySpan, zSpan, margin)