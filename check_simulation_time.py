"""
check simulation and calculation time for treeseqs and sumstats calculation
"""


from os import listdir
from os.path import isfile, join
from datetime import datetime
import numpy as np


def list_files(mypath):
    """just to list files"""
    return [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]


def grep_duration_from_logfile(path, year=2022):
    """greps times and calculates diffs"""
    times = [None]
    del times[0]

    try:
        with open(path, "r", encoding="utf-8") as infile:
            for line in infile:
                if line.startswith(str(year)):
                    line = line.rstrip()
                    times.append(line.split("\t")[0])

        if len(times) > 0:
            start = datetime.strptime(times[0], "%Y-%m-%d %H:%M:%S.%f")
            end = datetime.strptime(times[-1], "%Y-%m-%d %H:%M:%S.%f")

            diff = end - start
        else:
            diff = None

    except:
        diff = None

    return diff


def get_median_duration_from_logs(log_dir, qantiles=[0.025, 0.25, 0.5, 0.75, 0.975]):
    """calculate quantiles"""
    simts_files = list_files(log_dir)

    simts_duration = [None]
    del simts_duration[0]
    for file in simts_files:
        duration = grep_duration_from_logfile(file)
        simts_duration.append(duration)

    simts_duration = np.array(
        [duration for duration in simts_duration if not duration is None]
    )

    if len(simts_duration) > 0:
        quantiles = np.quantile(simts_duration, qantiles)
    else:
        quantiles = np.array([None for _ in quantiles])

    return quantiles


my_dirs = [
    "logs/module01/simulate_treeseqs/",
    "logs/module01/alternative_model_simulate_treeseqs/",
    "logs/module01/calculate_sumstats/",
    "logs/module01/alternative_model_calculate_sumstats/",
    "logs/module01/calculate_masked_sumstats/",
    "logs/module01/alternative_model_calculate_masked_sumstats/",
]


q = [0.025, 0.25, 0.5, 0.75, 0.975]
print(f"quantiles: {q}")
for this_dir in my_dirs:
    _ = [print(d, end="\t") for d in get_median_duration_from_logs(this_dir)]
    print(f"\t{this_dir}")
