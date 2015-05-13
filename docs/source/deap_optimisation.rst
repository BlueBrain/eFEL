DEAP optimisation
=================

.. contents::

Introduction
------------
Using the eFEL, pyNeuron and the DEAP optimisation library one can very easily 
set up a genetic algorithm to fit parameters of a neuron model.

We propose this setup because it leverages the power of the Python language 
to load several software tools in a compact script. The DEAP 
(Distributed Evolutionary Algorithms in Python) allows you to easily switch
algorithms. Parallelising your evaluation function over cluster computers 
becomes a matter of only adding a couple of lines to your 
`code <http://deap.readthedocs.org/en/latest/tutorials/basic/part4.html>`_, 
thanks to `pyScoop <http://pyscoop.org>`_.

In this example we will assume you have installed 
`eFEL <https://github.com/BlueBrain/eFEL>`_, 
`pyNeuron <http://www.neuron.yale.edu/neuron/download/compile_linux#otheroptions>`_ 
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

Evaluation function
-------------------
We now have to use pyNeuron to define the evaluation function to be optimised.
The input arguments are the parameters::

    g_pas, e_pas

and the return values::

    abs(voltage_base - target_voltage1)
    abs(steady_state_voltage - target_voltage2)

This translates into the following file (let's call it 'deap_efel_eval1.py')::

    import neuron
    neuron.h.load_file('stdrun.hoc')

    import efel

    # pylint: disable=W0212


    def evaluate(individual, target_voltage1=-80, target_voltage2=-60):
        """                                                                          
        Evaluate a neuron model with parameters e_pas and g_pas, extracts            
        eFeatures from resulting traces and returns a tuple with                      
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

Setting up the algorithm
------------------------
Now that we have an evaluation function we just have to pass this to the DEAP 
optimisation library. DEAP allows you to easily set up a genetic algorithm
to optimise your evaluation function. Let us first import all the necessary 
components::

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

Next we define a number of constants that will be used as settings for DEAP 
later::

    # Population size
    POP_SIZE = 100
    # Number of offspring in every generation
    OFFSPRING_SIZE = 100

    # Number of generations
    NGEN = 300

    # The parent and offspring population size are set the same
    MU = OFFSPRING_SIZE
    LAMBDA = OFFSPRING_SIZE
    # Crossover probability 
    CXPB = 0.7
    # Mutation probability, should sum to one together with CXPB
    MUTPB = 0.3

    # Eta parameter of cx and mut operators
    ETA = 10.0


We have two parameters with the following bounds::

    # The size of the individual is 2 (parameters g_pas and e_pas)
    IND_SIZE = 2

    LOWER = [1e-8, -100.0]
    UPPER = [1e-4, -20.0]


As evolutionary algorithm we choose 
`NSGA2 <http://www.tik.ee.ethz.ch/pisa/selectors/nsga2/nsga2_documentation.txt>`_::

    SELECTOR = "NSGA2"


Let's create the DEAP individual and fitness. 
We set the weights of the fitness values to -1.0 so that the fitness function 
will be minimised instead of maximised::

    creator.create("Fitness", base.Fitness, weights=[-1.0] * 2)

The individual will just be a list (of two parameters)::

    creator.create("Individual", list, fitness=creator.Fitness)

We want to start with individuals for which the parameters are picked from a 
uniform random distribution. Let's create a function that returns such a 
random list based on the bounds and the dimensions of the problem::

    def uniform(lower_list, upper_list, dimensions):
        """Fill array """

        if hasattr(lower_list, '__iter__'):
            return [random.uniform(lower, upper) for lower, upper in
                    zip(lower_list, upper_list)]
        else:
            return [random.uniform(lower_list, upper_list)
                    for _ in range(dimensions)]

DEAP works with the concept of 'toolboxes'. The user defines genetic 
algorithm's individuals, operators, etc by registering them in a toolbox.

We first create the toolbox::

    toolbox = base.Toolbox()

Then we register the 'uniform' function we defined above::

    toolbox.register("uniformparams", uniform, LOWER, UPPER, IND_SIZE)

The three last parameters of this register call will be passed on to the
'uniform' function call

Now we can also register an individual::

    toolbox.register(
        "Individual",
        tools.initIterate,
        creator.Individual,
        toolbox.uniformparams)

And a population as list of individuals::

    toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

The function to evaluate we defined above. Assuming you saved that files as
'deap_efel_eval1.py', we can import it as a module, and register the function::

    import deap_efel_eval1
    toolbox.register("evaluate", deap_efel_eval1.evaluate)

For the mutation and crossover operator we use builtin operators that are
typically used with NSGA2::

    toolbox.register(
        "mate",
        deap.tools.cxSimulatedBinaryBounded,
        eta=ETA,
        low=LOWER,
        up=UPPER)
    toolbox.register("mutate", deap.tools.mutPolynomialBounded, eta=ETA,
                     low=LOWER, up=UPPER, indpb=0.1)

And then we specify the selector to be used::

    toolbox.register(
        "select",
        tools.selNSGA2)

We initialise the population with the size of the offspring::

    pop = toolbox.population(n=MU)


And register some statistics we want to print during the run of the algorithm::

    first_stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    second_stats = tools.Statistics(key=lambda ind: ind.fitness.values[1])
    stats = tools.MultiStatistics(obj1=first_stats, obj2=second_stats)
    stats.register("min", numpy.min, axis=0)

The only thing that is left now is to run the algorithm in 'main'::

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

For you convenience the full code is in a code block below. It should be saved
as 'deap_efel.py'.

Running the code
----------------
Assuming that the necessary dependencies are installed correctly the 
optimisation can then be run with::

    python deap_efel.py

The full code of 'deap_efel.py'::

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
