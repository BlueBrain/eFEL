from __future__ import annotations
from unittest.mock import patch
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
        self.features = ["single_burst_ratio", "ISIs", "irregularity_index",
                         "ISI_log_slope", "ISI_semilog_slope", "ISI_log_slope_skip",
                         "burst_ISI_indices"
                         ]
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

    def test_irregularity_index(self):
        assert (
            self.feature_values["irregularity_index"] == pytest.approx(0.0, abs=1e-9)
        )

    def test_ISI_log_slope(self):
        assert self.feature_values["ISI_log_slope"] == pytest.approx(0.0)

    def test_ISI_semilog_slope(self):
        assert self.feature_values["ISI_semilog_slope"] == pytest.approx(0.0)

    def test_ISI_log_slope_skip(self):
        assert self.feature_values["ISI_log_slope_skip"] == pytest.approx(0.0)

    def test_ISI_log_slope_skip_ValueError(self):
        with pytest.raises(ValueError):
            efel.set_double_setting("spike_skipf", 1.0)
            get_feature_values(
                [self.trace],
                ["ISI_log_slope_skip"], raise_warnings=False)[0]

    def test_burst_ISI_indices(self):
        assert self.feature_values["burst_ISI_indices"] is None


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


class TestBurstISIIndices:
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        efel.reset()
        self.trace = {  # bypassed values since ISI_values are defined in the mock
            "T": [1, 2, 3],
            "V": [4, 5, 6],
            "stim_start": [0],
            "stim_end": [1],
        }
        self.features = ["burst_ISI_indices"]
        yield
        efel.reset()

    @patch('efel.pyfeatures.isi.ISI_values')
    def test_burst_ISI_indices(self, mock_ISI_values):
        mock_ISI_values.return_value = np.array([50, 100, 50, 200, 50])

        self.feature_values = efel.getFeatureValues(
            [self.trace], self.features, raise_warnings=False
        )[0]

        assert self.feature_values["burst_ISI_indices"].tolist() == [2, 4]
