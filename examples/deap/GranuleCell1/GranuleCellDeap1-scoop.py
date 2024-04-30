"""Neuron / DEAP optimisation example using a cerebellar Granule cell model (D'Angelo et al. 2001)"""

# These scripts were developed (partly based on existing code) during the 2015 HBP Summit Hackathon
# by
# Werner Van Geit @ BlueBrain Project
# Stefan Masoli @ Department of Brain and Behavioral Sciences, Pavia University, Pavia

import random

import neuron
neuron.h.load_file('stdrun.hoc')

import efel
import numpy as np
import GranuleCell1 as grcsteps
import bluepyopt.deapext

from scoop import futures

# dictionary with all the original conductances
import collections
orig_conductances = collections.OrderedDict()
orig_conductances['GrC_CaHVA'] = 0.00046
orig_conductances['GrC_KA'] = 0.004
orig_conductances['GrC_Lkg1'] = 5.68e-5
orig_conductances['GrG_KM'] = 0.00035
orig_conductances['Calc'] = 1.5
orig_conductances['GrC_KCa'] = 0.004
orig_conductances['GrC_Kir'] = 0.0009
orig_conductances['GrC_Lkg2'] = 2.17e-5
orig_conductances['GrC_pNa'] = 2e-5
orig_conductances['GrG_KV'] = 0.003
orig_conductances['GrG_Na'] = 0.013
orig_conductances['GrG_Nar'] = 0.0005

import deap.gp
from deap import base
from deap import creator
from deap import tools

FEATURES = ['voltage_base',
            'spike_count',
            'min_voltage_between_spikes',
            'steady_state_voltage_stimend']

def get_features(conductances):
    feature_values = collections.defaultdict(dict)
    traces = []
    times = []
    voltages = []
    for step_number in range(3):
        time, voltage = grcsteps.steps(step_number, conductances)

        trace = {}  # dictionary with time, voltage start and end vectors.
        trace['T'] = time
        trace['V'] = voltage
        trace['stim_start'] = [100]
        trace['stim_end'] = [500]
        traces = [trace]  # a think is a list from a dictionary

        result = efel.get_feature_values(traces, FEATURES, raise_warnings=False)

        for feature_name, feature_list in result[0].items():
            if feature_list is not None and len(feature_list) > 0:
                feature_values[step_number][
                    feature_name] = np.mean(feature_list)
            else:
                feature_values[step_number][feature_name] = None

        times.append(time)
        voltages.append(voltage)

    return feature_values, times, voltages

orig_features, orig_times, orig_voltages = get_features(orig_conductances)

def evaluate(individual):
    """
    Evaluate a neuron model with parameters e_pas and g_pas, extracts
    features from resulting traces and returns a tuple with
    abs(voltage_base-target_voltage1) and
    abs(steady_state_voltage-target_voltage2)
    """

    conductances = dict()

    for conductance_number, key in enumerate(orig_conductances.keys()):
        conductances[key] = individual[conductance_number]

    new_features, times, voltages = get_features(conductances)
    individual.times = times
    individual.voltages = voltages
    obj_values = compare_features(orig_features, new_features)

    individual.tot_fitness = sum(obj_values)

    return obj_values


def plot_selector(population, k, selector=None, **kwargs):
    selected_pop = selector(population, k, **kwargs)
    return selected_pop

def compare_features(orig_features, new_features):
    obj_values = []
    for step_number in orig_features:
        for feature_name, feature_value in orig_features[step_number].items():
            if feature_value is not None:
                # print new_features
                # print step_number, feature_name
                if new_features[step_number][feature_name] is None:
                    obj_values.append(250.0)
                else:
                    obj_values.append(
                        abs(new_features
                            [step_number][feature_name] -
                            feature_value))
    return obj_values

random.seed(1)
POP_SIZE = 50
OFFSPRING_SIZE = 50
NGEN = 100
ALPHA = POP_SIZE
MU = OFFSPRING_SIZE
LAMBDA = OFFSPRING_SIZE

# Crossover, mutation probabilities
CXPB = 0.7

# Eta parameter of crossover / mutation parameters
# Basically defines how much they 'spread' solution around
# The higher this value, the more spread
ETA = 10.0

# Number of parameters
IND_SIZE = 12

# Bounds for the parameters
LOWER, UPPER = zip(*[(value - abs(value), value + abs(value))
                     for value in orig_conductances.values()])

# Number of objectives
OBJ_SIZE = 11

# Create a fitness function
# By default DEAP selector will try to optimise fitness values,
# so we add a -1 weight value to minise
creator.create("Fitness", base.Fitness, weights=[-1.0] * OBJ_SIZE)

class Individual(list):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.time = None
        self.voltage = None

# Create an individual that consists of a list
creator.create("Individual", Individual, fitness=creator.Fitness)

# Define a function that will uniformly pick an individual
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
toolbox.register("Individual", tools.initIterate, creator.Individual, toolbox.uniformparams)
toolbox.register("population", tools.initRepeat, list, toolbox.Individual)
toolbox.register("evaluate", evaluate)
toolbox.register("mate", deap.tools.cxSimulatedBinaryBounded, eta=ETA, low=LOWER, up=UPPER)
toolbox.register("mutate", deap.tools.mutPolynomialBounded, eta=ETA, low=LOWER, up=UPPER, indpb=0.1)
toolbox.register("variate", deap.algorithms.varAnd)
toolbox.register("select", plot_selector, selector=bluepyopt.deapext.tools.selIBEA)
toolbox.register("map", futures.map)

def main():
    pop = toolbox.population(n=MU)
    pop, logbook = bluepyopt.deapext.algorithms.eaAlphaMuPlusLambdaCheckpoint(pop, toolbox, MU,
                                                   CXPB, 1 - CXPB,
                                                  NGEN)
    return pop, logbook

if __name__ == '__main__':
    main()
