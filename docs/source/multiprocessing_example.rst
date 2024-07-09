Extracating Features in parallel using multiprocessing/scoop
============================================================

You can use multiprocessing or scoop (Scalable COncurrent Operations in
Python) packages to extract features in parallel. This is a simple
example of how to do it using multiprocessing:

.. code:: ipython3

    import efel
    import numpy

.. code:: ipython3

    def main():
        """Main"""
    
        traces = []
        for filename in ['example_trace1.txt', 'example_trace2.txt']:
            # Use numpy to read the trace data from the txt file
            data = numpy.loadtxt(filename)
    
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
            traces.append(trace1)
    
        # Now we pass 'traces' to the efel and ask it to calculate the feature
        # values
    
        ## Using mutliprocessing
        import multiprocessing
        pool = multiprocessing.Pool()
        traces_results = efel.get_feature_values(
            traces, [
                'AP_amplitude', 'voltage_base'], parallel_map=pool.map)
    
        # The return value is a list of trace_results, every trace_results
        # corresponds to one trace in the 'traces' list above (in same order)
        for trace_number, trace_results in enumerate(traces_results):
            print("Results for trace %d" % (trace_number + 1))
            # trace_result is a dictionary, with as keys the requested features
            for feature_name, feature_values in trace_results.items():
                print("Feature %s has the following values: %s" % \
                    (feature_name, ', '.join([str(x) for x in feature_values])))

.. code:: ipython3

    if __name__ == '__main__':
        main()

If you want to use scoop, it is advisable to use a python script
scoop\_example.py (present in the ``parallel`` example folder) instead
of the notebook e.g. using

``python -m scoop -n 4 scoop_example.py``

where arguments are:

``-m scoop``: The ``-m`` flag tells Python to run a module as a script.
In this case, the module is scoop, which is a library for distributed
computing. Running SCOOP as a module initializes its environment,
preparing it for distributed execution.

``-n 4``: This option specifies the number of worker processes that
scoop should use for executing tasks. In this case, 4 workers are
requested. scoop will distribute tasks among these workers, allowing for
parallel execution.
