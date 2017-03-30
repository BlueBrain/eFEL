
import neo
import efel
import quantities as pq


#first file
file_name = "neo_test_file_no_times.pickle"
bl = neo.core.Block()
seg = neo.core.Segment()
data = range(10)
rate = 1000*pq.Hz
signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
seg.analogsignals.append(signal)
bl.segments.append(seg)
neo.io.PickleIO(filename=file_name, mode="rw").write_block(bl)

#second with times in epoch 
file_name = "neo_test_file_epoch_times.pickle"
bl = neo.core.Block()
seg = neo.core.Segment()
data = range(10)
rate = 1000*pq.Hz
signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
seg.analogsignals.append(signal)
seg.epochs.append(neo.core.Epoch(times=pq.Quantity([0.0,20.0], units=pq.ms), name="stim"))
bl.segments.append(seg)
neo.io.PickleIO(filename=file_name, mode="rw").write_block(bl) 

#events complete
file_name = "neo_test_file_events_time.pickle"
bl = neo.core.Block()
seg = neo.core.Segment()
data = range(10)
rate = 1000*pq.Hz
signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
seg.analogsignals.append(signal)
seg.events.append(neo.core.Event(times=[0.0]*pq.ms, units=pq.ms, name="stim_start" ))
seg.events.append(neo.core.Event(times=[20.0]*pq.ms, units=pq.ms, name="stim_end" ))
bl.segments.append(seg)
neo.io.PickleIO(filename=file_name, mode="rw").write_block(bl)

#events incomplete
file_name = "neo_test_file_events_time_incomplete.pickle"
bl = neo.core.Block()
seg = neo.core.Segment()
data = range(10)
rate = 1000*pq.Hz
signal = neo.core.AnalogSignal(data, sampling_rate=rate, units="mV")
seg.analogsignals.append(signal)
seg.events.append(neo.core.Event(times=[0.0]*pq.ms, units=pq.ms , name="stim_start" ))
bl.segments.append(seg)
neo.io.PickleIO(filename=file_name, mode="rw").write_block(bl)



