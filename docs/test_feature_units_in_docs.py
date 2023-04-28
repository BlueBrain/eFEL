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
