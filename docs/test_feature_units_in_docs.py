"""Parse eFeatures.rst.
    Make sure every unit there is present in efel.units
"""

from pathlib import Path
import re


from efel import units


def extract_feature_and_units(rst_text: str) -> dict:
    """Parse rst to extract features and units."""
    # Use regex to match and extract the desired information
    pattern = r'LibV\d+\s*:\s*([\w_]+)[\s\S]*?- \*\*Units\*\*:\s*([\w/\(\)]+)'
    matches = re.findall(pattern, rst_text)

    result = {}
    if matches:
        for match in matches:
            feature_name, unit = match
            result[feature_name] = unit
    return result


def test_extract_feature_and_units():
    rst1 = """
        LibV5 : voltage_base
        ~~~~~~~~~~~~~~~~~~~~

        The average voltage during the last 10% of time before the stimulus.

        - **Required features**: t, V, stim_start, stim_end
        - **Parameters**: voltage_base_start_perc (default = 0.9), voltage_base_end_perc (default = 1.0)
        - **Units**: mV
        - **Pseudocode**: ::

            voltage_base = numpy.mean(voltage[numpy.where(
                (t >= voltage_base_start_perc * stim_start) &
                (t <= voltage_base_end_perc * stim_start))])
        """
    assert extract_feature_and_units(rst1) == {"voltage_base": "mV"}

    rst2 = """
        LibV2 : AP_fall_rate_change
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Difference of the fall rates of the second and the first action potential
        divided by the fall rate of the first action potential

        - **Required features**: LibV2: AP_fall_rate_change
        - **Units**: constant
        - **Pseudocode**: ::

            AP_fall_rate_change = (AP_fall_rate[1:] - AP_fall_rate[0]) / AP_fall_rate[0]

        LibV5 : AP_phaseslope
        ~~~~~~~~~~~~~~~~~~~~~

        Slope of the V, dVdt phasespace plot at the beginning of every spike

        (at the point where the derivative crosses the DerivativeThreshold)

        - **Required features**: LibV5:AP_begin_indices
        - **Parameters**: AP_phaseslope_range
        - **Units**: 1/(ms)
        - **Pseudocode**: ::

            range_max_idxs = AP_begin_indices + AP_phseslope_range
            ...
    """
    assert extract_feature_and_units(rst2) == {
        "AP_fall_rate_change": "constant",
        "AP_phaseslope": "1/(ms)",
    }
    rst3 = """lorem ipsum dolor sit amet, consectetur adipiscing elit."""
    assert extract_feature_and_units(rst3) == {}

    rst4 = """
        LibV5 : AP_peak_upstroke
        ~~~~~~~~~~~~~~~~~~~~~~~~

        Maximum of rise rate of spike

        - **Required features**: LibV5: AP_begin_indices, LibV5: peak_indices
        - **Units**: V/s"""
    assert extract_feature_and_units(rst4) == {"AP_peak_upstroke": "V/s"}


def test_efeature_units():
    """Test to assure the rst and efel API are consistent."""
    # Read the rst file
    rst_file = Path(__file__).parent / "source" / "eFeatures.rst"
    with open(rst_file, "r") as f:
        rst_text = f.read()

    # Extract the feature names and units
    feature_units = extract_feature_and_units(rst_text)

    # Compare with the units in efel.units
    for feature_name, unit in feature_units.items():
        assert units.get_unit(feature_name) == unit
