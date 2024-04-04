"""Contains the features that are computed using multiple traces."""


import numpy as np

import efel


def bpap_attenuation(soma_trace: dict, dendrite_trace: dict) -> float:
    """Computes the attenuation of backpropagating action potential.

    Backpropagating action potential is the action potential that is initiated
    in the soma and propagates to the dendrite. The attenuation is the ratio
    of the amplitude of the action potential in the soma and the dendrite.
    The attenuation is computed by first subtracting the resting potential
    from the voltage traces.
    """
    f_values = efel.get_feature_values(
        [soma_trace, dendrite_trace], ["voltage_base"],
        parallel_map=None, return_list=True, raise_warnings=True)
    vb_soma = f_values[0]["voltage_base"][0]
    vb_dend = f_values[1]["voltage_base"][0]
    v_soma = soma_trace["V"]
    v_dend = dendrite_trace["V"]
    res = (np.max(v_soma) - vb_soma) / (np.max(v_dend) - vb_dend)
    return res
