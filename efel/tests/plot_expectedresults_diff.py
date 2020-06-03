#!/usr/bin/env python

"""Compare the current expected results file with the old expected
    results and plot the difference. """


"""
Copyright (c) 2015, Blue Brain Project/EPFL

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


import pandas as pd
import json
from scipy.spatial import distance
import matplotlib.pyplot as plt
from tqdm import tqdm
from pprint import pprint
import argparse

plt.style.use('ggplot')

parser = argparse.ArgumentParser(
    description='Plot the difference in expected results.')
parser.add_argument('old_path', type=str,
                    help='path to the old expected results json file')
args = parser.parse_args()

with open(args.old_path) as f:
    old_expected_results = json.load(f)

with open('testdata/allfeatures/expectedresults.json') as f:
    new_expected_results = json.load(f)


assert(set(old_expected_results.keys()) == set(new_expected_results.keys()))


diff_results = {}
for key, old_value in tqdm(old_expected_results.items()):
    if old_value in [None, []]:
        if new_expected_results[key] not in [None, []]:
            print(f"{key} was None or [] previously, now it has a value")
        difference = old_value
    else:
        difference = distance.euclidean(new_expected_results[key], old_value)
    diff_results[key] = difference


changed_results = {k: v for k, v in diff_results.items() if v and v > 0}


print("The difference is:")
pprint(changed_results)

df_diff = pd.DataFrame(changed_results, index=["distance"])


ax = df_diff.plot.bar(figsize=(16, 8),
                      title="Difference in expected results (Euclidean)")
for i in ax.patches:
    ax.text(i.get_x() + .025, i.get_height() + .0002,
            str(round((i.get_height()), 4)), color='dimgrey', rotation=0)
plt.savefig("diff_expected_results.png")
