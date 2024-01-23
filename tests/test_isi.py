import numpy as np
import pytest
import efel
from efel import get_feature_names, get_feature_values


def generate_spike_data(
    duration=1000, spike_interval=100, spike_duration=1, baseline=-80, spike_value=0
) -> tuple[np.ndarray, np.ndarray]:
    """Generate time and voltage arrays for a signal with spikes."""
    voltage = np.full(duration, baseline)
    # Add spikes to the voltage array
    for i in range(0, duration, spike_interval):
        voltage[i: i + spike_duration] = spike_value
    time = np.arange(duration)
    return time, voltage


class TestRegularISI:
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        # setup
        efel.reset()
        self.time, self.voltage = generate_spike_data()
        self.trace = {}
        self.trace['T'] = self.time
        self.trace['V'] = self.voltage
        self.trace['stim_start'] = [0]
        self.trace['stim_end'] = [1000]
        self.features = ["single_burst_ratio", "ISIs"]
        self.feature_values = get_feature_values(
            [self.trace],
            self.features, raise_warnings=False)[0]
        # teardown
        yield
        efel.reset()

    def test_ISIs(self):
        assert np.allclose(self.feature_values["ISIs"], 100.0)

    def test_single_burst_ratio(self):
        assert "single_burst_ratio" in get_feature_names()
        assert self.feature_values["single_burst_ratio"] == pytest.approx(1.0)


class TestThreeSpikes:
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        # setup
        efel.reset()
        duration = 100
        spike_interval = 25
        self.time, self.voltage = generate_spike_data(duration, spike_interval)
        self.trace = {}
        self.trace['T'] = self.time
        self.trace['V'] = self.voltage
        self.trace['stim_start'] = [0]
        self.trace['stim_end'] = [duration]
        self.features = ["single_burst_ratio", "ISIs", "peak_time"]
        self.feature_values = get_feature_values(
            [self.trace],
            self.features, raise_warnings=False)[0]
        # teardown
        yield
        efel.reset()

    def test_ISIs(self):
        # duration/spike_interval is 4 but the beginning of the voltage
        # trace is not a spike
        assert len(self.feature_values["peak_time"]) == 3
        assert len(self.feature_values["ISIs"]) == 2
        assert np.allclose(self.feature_values["ISIs"], 25.0)

    def test_single_burst_ratio(self):
        # none due to default ignore_first_ISI=True
        assert self.feature_values["single_burst_ratio"] is None

        # set ignore_first_ISI=False
        efel.set_int_setting("ignore_first_ISI", 0)
        self.feature_values = get_feature_values(
            [self.trace],
            self.features, raise_warnings=False)[0]
        assert self.feature_values["single_burst_ratio"] == pytest.approx(1.0)
