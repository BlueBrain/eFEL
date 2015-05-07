import efel
import json

efel.cppcore.Initialize(efel.settings.dependencyfile_path, "log")

with open('featurenames.json', 'w') as featurenames_json:
    feature_names = []
    efel.cppcore.getFeatureNames(feature_names)
    json.dump(
        feature_names,
        featurenames_json,
        indent=4,
        separators=(
            ',',
            ': '))
