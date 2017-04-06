Reading different file formats
==============================

Neo is a Python package which provides support for reading a wide range of neurophysiology file
formats, including Spike2, NeuroExplorer, AlphaOmega, Axon, Blackrock, Plexon and Tdt.

The function :func:`efel.io.load_neo_file()` reads data from any of the file formats supported by
Neo and formats it for use in eFEL.

As an example, suppose we have an .abf file containing a single trace. Since eFEL requires
information about the start and end times of the current injection stimulus, we provide these
times as well as the filename::

   import efel
   
   data = efel.io.load_neo_file("path/first_file.abf", stim_start=200, stim_end=700)

Since some file formats can contain multiple recording episodes (e.g. trials) and multiple
signals per episode, the function returns traces in a list of lists, like this::

   data : [Segment_1, Segment_2, ..., Segment_n]
            with Segment_1 = [Trace_1, Trace_2, ..., Trace_n]

Since our file contains only a single recording episode, our list of traces is::

   traces = data[0]

which we pass to eFEL as follows::

   features = efel.getFeatureValues(traces, ['AP_amplitude', 'voltage_base'])

Stimulus information within the file
------------------------------------

Some file formats can store information about the current injection stimulus. In this second
example, the file contains an :class:`Epoch` object named "stimulation", so we don't need to
explicitly specify `stim_start` and `stim_end`::

   data2 = efel.io.load_neo_file("path/second_file.h5")
