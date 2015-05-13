Quick start
===========

First you need to import the module::
    
    import efel
                                                                                 
To get a list with all the available eFeature names::

    efel.getFeatureNames()

The python function to extract eFeatures is getFeatureValues(...).
Below is a short example on how to use this function.

The code and example trace are available 
`here <https://github.com/BlueBrain/eFEL/blob/master/examples/basic/basic_example1.py>`_::

    """Basic example 1 for eFEL"""                                                   
                                                                                     
    import efel                                                                      
    import numpy                                                                     
                                                                                     
    def main():                                                                      
        """Main"""                                                                   
                                                                                     
        # Use numpy to read the trace data from the txt file                         
        data = numpy.loadtxt('example_trace1.txt')                                   
                                                                                     
        # Time is the first column                                                   
        time = data[:, 0]                                                            
        # Voltage is the second column                                               
        voltage = data[:, 1]                                                         
                                                                                     
        # Now we will construct the datastructure that will be passed to eFEL        
                                                                                     
        # A 'trace' is a dictionary                                                  
        trace1 = {}                                                                  
                                                                                     
        # Set the 'T' (=time) key of the trace                                       
        trace1['T'] = time                                                           
                                                                                     
        # Set the 'V' (=voltage) key of the trace                                    
        trace1['V'] = voltage                                                        
                                                                                     
        # Set the 'stim_start' (time at which a stimulus starts, in ms)              
        # key of the trace                                                           
        # Warning: this need to be a list (with one element)                         
        trace1['stim_start'] = [700]                                                 
                                                                                     
        # Set the 'stim_end' (time at which a stimulus end) key of the trace         
        # Warning: this need to be a list (with one element)                         
        trace1['stim_end'] = [2700]                                                  
                                                                                     
        # Multiple traces can be passed to the eFEL at the same time, so the         
        # argument should be a list                                                  
        traces = [trace1]                                                            
                                                                                     
        # Now we pass 'traces' to the efel and ask it to calculate the feature       
        # values                                                                     
        traces_results = efel.getFeatureValues(traces,                               
                                               ['AP_amplitude', 'voltage_base'])     
                                                                                     
        # The return value is a list of trace_results, every trace_results           
        # corresponds to one trace in the 'traces' list above (in same order)        
        for trace_results in traces_results:                                         
            # trace_result is a dictionary, with as keys the requested eFeatures      
            for feature_name, feature_values in trace_results.items():              
                print "Feature %s has the following values: %s" % \                  
                    (feature_name, ', '.join([str(x) for x in feature_values]))      
                                                                                     
                                                                                     
    if __name__ == '__main__':                                                       
        main()

The output of this example is::

    Feature AP_amplitude has the following values: 72.5782441262, 46.3672552618, 41.1546679158, 39.7631750953, 36.1614653031, 37.8489295737
    Feature voltage_base has the following values: -75.446665721

This means that the eFEL found 5 action potentials in the voltage trace. The     
amplitudes of these APs are the result of the 'AP_amplitude' feature.

The voltage before the start of the stimulis is measured by 'voltage_base'.      

Results are in mV.
