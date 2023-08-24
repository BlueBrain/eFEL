import json
import collections

er = collections.OrderedDict()

with open('expectedresults.txt', 'r') as er_file:
    er_content = er_file.read()

for line in er_content.split('\n'):
    split_line = line.split()
    if len(split_line) == 0:
        continue
    elif len(split_line) == 1:
        er[split_line[0]] = None
    else:
        er[split_line[0]] = [str(x) for x in split_line[1:]]

with open('expectedresults.json', 'w') as er_json:
    json.dump(er, er_json, indent=4, separators=(',', ': '))
