"""eFEL / DEAP / Neuron optimisation example"""

import random
import numpy

import deap
import deap.gp
import deap.benchmarks
from deap import base
from deap import creator
from deap import tools
from deap import algorithms

# Set random seed
random.seed(1)

# Set number of individuals in population
POP_SIZE = 100

# Set number of individuals to create in every offspring
OFFSPRING_SIZE = 100

# Total number of generation to run
NGEN = 300

# Total population size of EA
ALPHA = POP_SIZE
# Total parent population size of EA
MU = OFFSPRING_SIZE
# Total offspring size of EA
LAMBDA = OFFSPRING_SIZE

# Crossover, mutation probabilities
CXPB = 0.7
MUTPB = 0.3

# Eta parameter of crossover / mutation parameters
# Basically defines how much they 'spread' solution around
# The higher this value, the more spread
ETA = 10.0

# Selector to use
SELECTOR = "NSGA2"

# Number of parameters
IND_SIZE = 2

# Bounds for the parameters
LOWER = [1e-8, -100.0]
UPPER = [1e-4, -20.0]

# Create a fitness function
# By default DEAP selector will try to optimise fitness values,
# so we add a -1 weight value to minimise
creator.create("Fitness", base.Fitness, weights=[-1.0] * 2)

# Create an individual that consists of a list
creator.create("Individual", list, fitness=creator.Fitness)


# Define a function that will uniformly pick an individual
def uniform(lower_list, upper_list, dimensions):
    """Fill array """

    if hasattr(lower_list, '__iter__'):
        return [random.uniform(lower, upper) for lower, upper in
                zip(lower_list, upper_list)]
    else:
        return [random.uniform(lower_list, upper_list)
                for _ in range(dimensions)]

# Create a DEAP toolbox
toolbox = base.Toolbox()

# Register the 'uniform' function
toolbox.register("uniformparams", uniform, LOWER, UPPER, IND_SIZE)

# Register the individual format
# An indiviual is create by 'creator.Individual' and parameters are initially
# picked by 'uniform'
toolbox.register(
    "Individual",
    tools.initIterate,
    creator.Individual,
    toolbox.uniformparams)

# Register the population format. It is a list of individuals
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)

# Register the evaluation function for the individuals
import deap_efel_eval1
toolbox.register("evaluate", deap_efel_eval1.evaluate)

# Register the mate operator
toolbox.register(
    "mate",
    deap.tools.cxSimulatedBinaryBounded,
    eta=ETA,
    low=LOWER,
    up=UPPER)

# Register the mutation operator
toolbox.register("mutate", deap.tools.mutPolynomialBounded, eta=ETA,
                 low=LOWER, up=UPPER, indpb=0.1)

# Register the variate operator
toolbox.register("variate", deap.algorithms.varAnd)

# Register the selector (picks parents from population)
toolbox.register(
    "select",
    tools.selNSGA2)

# Generate the population object
pop = toolbox.population(n=MU)

# Register the statistics we want to record during the optimisation
# In this case only the minimal value
first_stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
second_stats = tools.Statistics(key=lambda ind: ind.fitness.values[1])
stats = tools.MultiStatistics(obj1=first_stats, obj2=second_stats)
stats.register("min", numpy.min, axis=0)

# Run the actual algorithm with all the object/parameters we set up above
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
