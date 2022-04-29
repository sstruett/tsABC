"""
check simulation and calculation time for treeseqs and sumstats calculation
"""

import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime

from os import listdir
from os.path import isfile, join


mypath = "logs/module01/simulate_treeseqs/"


def list_files(mypath):
    return [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]


def grep_duration_from_logfile(path, year=2022):
    times = []
    with open(path, "r") as infile:
        for line in infile:
            if line.startswith(str(year)):
                line = line.rstrip()
                times.append(line.split("\t")[0])

    start = datetime.strptime(times[0], "%Y-%m-%d %H:%M:%S.%f")
    end = datetime.strptime(times[-1], "%Y-%m-%d %H:%M:%S.%f")

    return end - start


def get_median_duration_from_logs(log_dir, q=[0.025, 0.25, 0.5, 0.75, 0.975]):
    simts_files = list_files(log_dir)

    simts_duration = []
    for file in simts_files:
        duration = grep_duration_from_logfile(file)
        simts_duration.append(duration)

    simts_duration = np.array(simts_duration)
    return np.quantile(simts_duration, q)


my_dirs = [
    "logs/module01/simulate_treeseqs/",
    "logs/module01/alternative_model_simulate_treeseqs",

    "logs/module01/calculate_sumstats/",
    "logs/module01/alternative_model_calculate_sumstats",
    
    "logs/module01/calculate_masked_sumstats",
    "logs/module01/alternative_model_calculate_masked_sumstats",
]


q=[0.025, 0.25, 0.5, 0.75, 0.975]
print(f"quantiles: {q}")
for this_dir in my_dirs:
    _ = [print(d, end="\t") for d in get_median_duration_from_logs(this_dir)]
    print(f"\t{this_dir}")





