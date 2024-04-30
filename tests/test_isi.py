from __future__ import annotations
from unittest.mock import MagicMock, patch
import numpy as np
import pytest
import efel
from efel import get_feature_names, get_feature_values
from efel.pyfeatures import isi


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


def test_ISIs_single_spike():
    """Test the edge case where there are less than 2 spikes."""
    time, voltage = generate_spike_data(spike_interval=600)
    trace = {
        "T": time,
        "V": voltage,
        "stim_start": [0],
        "stim_end": [1000],
    }

    features = ["peak_time", "ISIs"]
    feature_values = efel.get_feature_values([trace], features)[0]

    assert len(feature_values["peak_time"]) == 1
    assert feature_values["ISIs"] is None


def test_isi_log_slope_core_exception():
    """Test that _isi_log_slope_core handles np.polyfit raising a LinAlgError."""
    # Mock np.polyfit to raise a LinAlgError
    np.polyfit = MagicMock(side_effect=np.linalg.LinAlgError("Singular matrix"))

    # Call _isi_log_slope_core with some dummy data
    isi_values = np.array([1, 2, 3, 4, 5])
    with pytest.warns(UserWarning) as record:
        result = isi._isi_log_slope_core(isi_values)

    assert "Error in polyfit: Singular matrix" in str(record[0].message)
    assert result is None

    # Reset np.polyfit to its original state
    np.polyfit = np.lib.polynomial.polyfit


class TestRegularISI:
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        # setup
        efel.reset()
        self.time, self.voltage = generate_spike_data()
        self.trace = {
            'T': self.time,
            'V': self.voltage,
            'stim_start': [0],
            'stim_end': [1000],
        }
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
        self.trace = {
            'T': self.time,
            'V': self.voltage,
            'stim_start': [0],
            'stim_end': [duration],
        }
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

        self.feature_values = efel.get_feature_values(
            [self.trace], self.features, raise_warnings=False
        )[0]

        assert self.feature_values["burst_ISI_indices"].tolist() == [2, 4]
