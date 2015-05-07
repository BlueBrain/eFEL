"""efel main python functions"""

"""
Copyright (c) 2015, EPFL/Blue Brain Project

 This file is part of eFEL <https://github.com/BlueBrain/eFEL>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

# pylint: disable=W0602,W0603,W0702, F0401, W0612, R0912

import sys
sys.path.append("/home/vangeit/local/Feature/lib")
import felipy
import numpy
import os
import fel

settings = fel.Settings()


def printelements(iterat):
    """Print all the elements of an iterator"""
    return str([element for element in iterat])


def setDependencyFileLocation(location):
    """Set the location of the DependencyV...txt file"""

    global dependencyFileLocation
    if not os.path.exists(location):
        raise Exception(
            "Path to dependency file {%s} doesn't exist" %
            location)
    settings.dependencyfile_path = location


def getDependencyFileLocation():
    """Get the location of the DependencyV...txt file"""

    return settings.dependencyfile_path


def setThreshold(newThreshold):
    """Set the spike detection threshold in the FEL"""
    settings.threshold = newThreshold


def setDerivativeThreshold(newDerivativeThreshold):
    """Set the threshold for the derivate for detecting the spike onset"""
    settings.derivate_threshold = newDerivativeThreshold


def getFeatureValues(traces, featureNames):
    """
    input :
       traces: dictionary with keys 'T' (time values list),
               'V' (voltage values list), 'stim_start', stim_end'
       featureNames: list of strings
    output :
       featureDicts: list of dictionaries (1 for each trace),
                     featureNames are the keys in the dicts, the values in the
                     dicts are lists with detected feature values per trace
    """
    featureDicts = []

    for trace in traces:
        featureDict = {}

        if 'stim_start' in trace and 'stim_end' in trace:
            if len(trace['stim_start']) == 1 and len(trace['stim_end']) == 1:
                if trace['stim_end'][0] <= trace['stim_start'][0]:
                    raise Exception(
                        'stim_end needs to be larger than '
                        'stim_start:\nstim_start=%f stim_end=%f' %
                        (trace['stim_start'][0], trace['stim_end'][0]))
            else:
                raise Exception('stim_start and stim_end in the trace '
                                'dictionary need to be lists of 1 element')

        else:
            raise Exception('stim_start or stim_end missing from trace')

        felipy.Initialize(settings.dependencyfile_path, "log")

        # First set some settings that are used by the feature extraction
        felipy.setFeatureDouble('spike_skipf', [0.1])
        felipy.setFeatureInt('max_spike_skip', [2])
        felipy.setFeatureDouble('Threshold',
                                [settings.threshold])
        felipy.setFeatureDouble('DerivativeThreshold',
                                [settings.derivative_threshold])
        felipy.setFeatureDouble('interp_step', [0.1])
        felipy.setFeatureDouble('burst_factor', [1.5])
        felipy.setFeatureDouble("initial_perc", [0.1])

        # Next set time, voltage and the stimulus start and end
        for item in trace.keys():
            felipy.setFeatureDouble(item, [x for x in trace[item]])

        for featureName in featureNames:
            featureType = felipy.featuretype(featureName)
            if featureType == "double":
                felipyFeatureValues = list()
                try:
                    exitCode = felipy.getFeatureDouble(
                        featureName,
                        felipyFeatureValues)
                except:
                    exitCode = -1
            elif featureType == "int":
                felipyFeatureValues = list()
                exitCode = felipy.getFeatureInt(
                    featureName,
                    felipyFeatureValues)
            else:
                print "Feature %s has an unknown type: %s" % \
                    (featureName, featureType)
                exit(1)
            if exitCode < 0:
                import warnings
                warnings.warn(
                    "Error while calculating feature %s: %s" %
                    (featureName, felipy.getgError()),
                    RuntimeWarning)
                featureDict[featureName] = None
            else:
                featureDict[featureName] = numpy.array(felipyFeatureValues)

        featureDicts.append(featureDict)
    return featureDicts


def getMeanFeatureValues(traces, featureNames):
    """
    input :
       traces: dictionary with keys 'T' (time values list),
               'V' (voltage values list), 'stim_start', stim_end'
       featureNames: list of strings
    output :
       featureDicts: list of dictionaries (1 for each trace),
                     featureNames are the keys in the dicts,
                     the values in the dicts are mean values per feature
    """
    featureDicts = getFeatureValues(traces, featureNames)
    for featureDict in featureDicts:
        for (key, values) in featureDict.items():
            featureDict[key] = numpy.mean(values)

    return featureDicts


def main():
    """Test main() function"""
    setDependencyFileLocation("DependencyV5.txt")
    featureNames = [
        "burst_number",
        "mean_frequency",
        "voltage_base",
        "AP_height",
        "time_to_first_spike",
        "adaptation_index",
        "spike_half_width",
        "AHP_depth_abs",
        "AHP_depth_diff"]

    fileNames = ["../test/allfeaturetest/testdata.txt"]

    traces = []
    for fileName in fileNames:
        trace = {}
        loadedTrace = numpy.loadtxt(fileName)
        trace['T'] = loadedTrace[:, 0]
        trace['V'] = loadedTrace[:, 1]
        print trace['T']
        trace['stim_start'] = [700]
        trace['stim_end'] = [2700]
        traces.append(trace)

    print getFeatureValues(traces, featureNames)

# This function is just for information purposes
if __name__ == "__main__":
    main()
