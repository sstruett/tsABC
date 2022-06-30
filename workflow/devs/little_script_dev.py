import sys
import datetime
import warnings
import itertools
import pickle
import math
import numpy as np
import pandas as pd
import tskit
import json
import pyfuncs  # from file


# read sample names
sample_names = np.loadtxt(snakemake.input.sample_list, dtype=str)


# read tree sequence
treeseq_athal = tskit.load(snakemake.input.athal_treeseq)


# find the node ids for the sample of the population
population_sample = []
for individual in treeseq_athal.individuals():
    if str(json.loads(individual.metadata)["id"]) in sample_names:
        population_sample.extend(individual.nodes)


# sample treeseq to provided samples
treeseq_athal_population = treeseq_athal.simplify(samples=population_sample)
del treeseq_athal


# get the chromosomal regions from the config file
chrom_multiplier = 1e8
chromosome_regions = []
for chromid, (start, stop) in enumerate(
    [
    (0, 30_427_671),
    (0, 19_698_289),
    (0, 23_459_830),
    (0, 18_585_056),
    (0, 26_975_451)
    ],
    start=1,
):
    start += chrom_multiplier * chromid
    stop += chrom_multiplier * chromid
    chromosome_regions.append((start, stop))


treeseq_list = []
for chromid, (start, stop) in enumerate(chromosome_regions):
    treeseq_list.append(treeseq_athal_population.keep_intervals([(start, stop)]).trim())


# calculate the diversity along the sequence
div_list = []
step = 100_000
for tsid, ts in enumerate(treeseq_list):
    windows = [0, ts.sequence_length]
    for i in range(int(ts.sequence_length)):
        if (i * step) > ts.sequence_length:
            break
        else:
            windows.append(i * step)
    windows = sorted(list(set(windows)))

    my_div = ts.diversity(windows=windows)

    div_list.append(list(zip(windows[0:-1], my_div)))
    print(f"done {tsid}")

df_list_raw = [pd.DataFrame(div) for div in div_list]

df_list = []
for dfid, df in enumerate(df_list_raw):
    df["chr"] = dfid + 1
    df_list.append(df)


df = pd.concat(df_list)

df.to_csv("diversity_sliding_window.csv")


