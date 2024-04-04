"""Neuron / DEAP optimisation example evaluation function"""

import neuron
neuron.h.load_file('stdrun.hoc')

import efel

# pylint: disable=W0212


def evaluate(individual, target_voltage1=-80, target_voltage2=-60):
    """
    Evaluate a neuron model with parameters e_pas and g_pas, extracts
    features from resulting traces and returns a tuple with
    abs(voltage_base-target_voltage1) and 
    abs(steady_state_voltage-target_voltage2)
    """

    neuron.h.v_init = target_voltage1

    soma = neuron.h.Section()

    soma.insert('pas')

    soma.g_pas = individual[0]
    soma.e_pas = individual[1]

    clamp = neuron.h.IClamp(0.5, sec=soma)

    stim_start = 500
    stim_end = 1000

    clamp.amp = 1.0
    clamp.delay = stim_start
    clamp.dur = 100000

    voltage = neuron.h.Vector()
    voltage.record(soma(0.5)._ref_v)

    time = neuron.h.Vector()
    time.record(neuron.h._ref_t)

    neuron.h.tstop = stim_end
    neuron.h.run()

    trace = {}
    trace['T'] = time
    trace['V'] = voltage
    trace['stim_start'] = [stim_start]
    trace['stim_end'] = [stim_end]
    traces = [trace]

    features = efel.get_feature_values(traces, ["voltage_base",
                                             "steady_state_voltage"])
    voltage_base = features[0]["voltage_base"][0]
    steady_state_voltage = features[0]["steady_state_voltage"][0]

    return abs(target_voltage1 - voltage_base), \
        abs(target_voltage2 - steady_state_voltage)
