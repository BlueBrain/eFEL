DEAP optimisation
=================

Using the eFEL, pyNeuron and the DEAP optimisation library one can very easily 
set up a genetic algorithm to fit parameters of a neuron model.

In this example we will assume you have installed 
`eFEL <https://github.com/BlueBrain/eFEL>`_, 
`pyNeuron <http://www.neuron.yale.edu/neuron/static/new_doc/index.html>`_ 
and `DEAP <https://github.com/DEAP/deap>`_

The code of the example below can be downloaded from 
`here <https://github.com/BlueBrain/eFEL/tree/master/examples/deap>`_

To keep the example simple, let's start from a passive single compartmental 
model. The parameters to fit will be the conductance and reversal potential 
of the leak channel. We will simulate the model for 1000 ms, and at 500 ms
a step current of 1.0 nA is injected until the end of the simulation.

The objective values of the optimisation will be the voltage before the 
current injection (i.e. the 'voltage_base' feature), and the steady state 
voltage during the current injection at the end of the simulation 
('steady_state_voltage').

We now have to use pyNeuron to define the evaluation function to be optimised.
The input arguments are the parameters::

    g_pas, e_pas

and the return values::

    abs(voltage_base - target_voltage1)
    abs(steady_state_voltage - target_voltage2)

This translates into::

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

        features = efel.getFeatureValues(traces, ["voltage_base",
                                                 "steady_state_voltage"])
        voltage_base = features[0]["voltage_base"][0]
        steady_state_voltage = features[0]["steady_state_voltage"][0]

        return abs(target_voltage1 - voltage_base), \
            abs(target_voltage2 - steady_state_voltage)

Now that we have an evaluation function we just have to pass this to DEAP::

    import random
    import numpy

    import deap
    import deap.gp
    import deap.benchmarks
    from deap import base
    from deap import creator
    from deap import tools
    from deap import algorithms

    random.seed(1)
    POP_SIZE = 100
    OFFSPRING_SIZE = 100

    NGEN = 300
    ALPHA = POP_SIZE
    MU = OFFSPRING_SIZE
    LAMBDA = OFFSPRING_SIZE
    CXPB = 0.7
    MUTPB = 0.3
    ETA = 10.0

    SELECTOR = "NSGA2"

    IND_SIZE = 2
    LOWER = [1e-8, -100.0]
    UPPER = [1e-4, -20.0]

    creator.create("Fitness", base.Fitness, weights=[-1.0] * 2)
    creator.create("Individual", list, fitness=creator.Fitness)


    def uniform(lower_list, upper_list, dimensions):
        """Fill array """

        if hasattr(lower_list, '__iter__'):
            return [random.uniform(lower, upper) for lower, upper in
                    zip(lower_list, upper_list)]
        else:
            return [random.uniform(lower_list, upper_list)
                    for _ in range(dimensions)]

    toolbox = base.Toolbox()
    toolbox.register("uniformparams", uniform, LOWER, UPPER, IND_SIZE)
    toolbox.register(
        "Individual",
        tools.initIterate,
        creator.Individual,
        toolbox.uniformparams)
    toolbox.register("population", tools.initRepeat, list, toolbox.Individual)


    import deap_efel_eval1
    toolbox.register("evaluate", deap_efel_eval1.evaluate)

    toolbox.register(
        "mate",
        deap.tools.cxSimulatedBinaryBounded,
        eta=ETA,
        low=LOWER,
        up=UPPER)
    toolbox.register("mutate", deap.tools.mutPolynomialBounded, eta=ETA,
                     low=LOWER, up=UPPER, indpb=0.1)

    toolbox.register("variate", deap.algorithms.varAnd)

    toolbox.register(
        "select",
        tools.selNSGA2)

    pop = toolbox.population(n=MU)

    first_stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    second_stats = tools.Statistics(key=lambda ind: ind.fitness.values[1])
    stats = tools.MultiStatistics(obj1=first_stats, obj2=second_stats)
    stats.register("min", numpy.min, axis=0)

    if __name__ == '__main__':
        pop, logbook = algorithms.eaMuPlusLambda(
            pop,
            toolbox,
            MU,
            LAMBDA,
            CXPB,
            MUTPB,
            NGEN,
            stats,
            halloffame=None)
