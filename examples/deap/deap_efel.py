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
