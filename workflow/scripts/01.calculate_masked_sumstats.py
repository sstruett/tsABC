"""
caluclate masked summarizing statistics on the simulated treesequences
"""


import sys
import datetime
import warnings
import itertools
import pickle
import math
import numpy as np
import tskit
import pyfuncs  # from file


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("start calculating summarizing statisticss", file=logfile)


# Loading list of tree sequences
with open(snakemake.input.tsl, "rb") as tsl_file:
    tsl = pickle.load(tsl_file)

# sample only one haplotype per individual
rng = np.random.default_rng(snakemake.params.seed)
tsl_haploid = []
for treeseq in tsl:
    sample_set = [
        this_individual.nodes[rng.integers(low=0, high=2)]
        for this_individual in treeseq.individuals()
    ]
    tsl_haploid.append(treeseq.simplify(samples=sample_set))
del tsl
tsl = tsl_haploid
del tsl_haploid


# read and prepare mask files
mask = np.loadtxt(snakemake.input.mask, dtype=int)
region_start, region_end = snakemake.config["ABC"]["athaliana"]["observations"][
    "treeseq_1001"
]["chosen_region"][int(float(snakemake.wildcards.locid))]
region_mask = mask[(region_start <= mask[:, 0]) & (region_end > mask[:, 1])]
mask = region_mask - region_start
del region_mask

# filter mask for disjoint intervals
mask = pyfuncs.filter_mask_for_disjoint_intervals(mask, log=snakemake.log.log1)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "proportion of exons to mask in chromosome(1-indexed) "
        + f"{int(float(snakemake.wildcards.locid)) + 1}: "
        + f"{(mask[:, 1] - mask[:, 0]).sum()/(region_end-region_start)}",
        file=logfile,
    )


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("read and prepared the mask file", file=logfile)


# mask for regions provided
tsl_masked = []
for treeseq in tsl:
    tsl_masked.append(
        treeseq.delete_intervals(mask, simplify=True, record_provenance=True)
    )
tsl = tsl_masked
del tsl_masked


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("read and masked the tree sequences", file=logfile)


# read the breakpoint files
for breaks_filename in snakemake.input.breakpoints:
    if "SFS" in breaks_filename:
        with warnings.catch_warnings():  # will prevent warning if file is empty
            warnings.simplefilter("ignore")
            breaks_sfs = np.loadtxt(breaks_filename, dtype=float)
    elif "LD" in breaks_filename:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            breaks_ld = np.loadtxt(breaks_filename, dtype=float)
    elif "TM_WIN" in breaks_filename:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            breaks_tm_win = np.loadtxt(breaks_filename, dtype=float)

# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("read the breaks for binning", file=logfile)


# which sumstats to calculate
listed_sumstats = set(
    itertools.chain.from_iterable(
        [sumstat_set.split("/") for sumstat_set in snakemake.params.sumstats]
    )
)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "calculating summarizing statistics",
        file=logfile,
    )

# calculate the SFS
if "SFS" in listed_sumstats:
    breaks = breaks_sfs
    dataframe_sumstats_sfs = np.array(
        [
            treeseq.allele_frequency_spectrum(
                sample_sets=None,
                windows=None,
                mode="site",
                span_normalise=False,
                polarised=snakemake.config["ABC"]["sumstats_specs"]["SFS"]["polarised"],
            )
            for treeseq in tsl
        ]
    )

    # check if the sfs shall be discretized as well; this usually is not
    # necessary as the sfs is a discrete statistic already
    if len(breaks):
        # for each sfs sum up the values in between the bin_edges (=breaks)
        assert (
            len(breaks) >= 2
        ), "expect closed interval for bin_edges; at least 2 breaks must be provided"
        new_sfs_dataframe = []
        for sfs in dataframe_sumstats_sfs:
            sfs_vals = []
            for breakid in range(1, len(breaks)):
                low = int(breaks[breakid - 1])
                high = breaks[breakid]
                if np.isinf(high):
                    high = len(sfs) + 1
                else:
                    high = math.ceil(high)
                sfs_vals.append(sum(sfs[low:high]))
            new_sfs_dataframe.append(sfs_vals)

        # reassign to the dataframe sumstat
        dataframe_sumstats_sfs = np.array(new_sfs_dataframe)
else:
    dataframe_sumstats_sfs = np.empty(shape=0)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "calculated SFS statistics",
        file=logfile,
    )


# calculate the LD
if "LD" in listed_sumstats:
    specs = snakemake.config["ABC"]["sumstats_specs"]["LD"]
    breaks = breaks_ld
    dataframe_sumstats_ld = np.array(
        # the rng is to sample down the number of sites if needed
        [
            pyfuncs.calculate_ld(treeseq, specs, breaks, rng, snakemake.log.log1)
            for treeseq in tsl
        ]
    )
else:
    dataframe_sumstats_ld = np.empty(shape=0)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "calculated LD statistics",
        file=logfile,
    )


# calculate the TM_WIN
if "TM_WIN" in listed_sumstats:
    specs = snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]
    breaks = breaks_tm_win
    dataframe_sumstats_tm_win = np.array(
        [pyfuncs.calculate_tm_win(treeseq, specs, breaks_tm_win) for treeseq in tsl]
    )
else:
    dataframe_sumstats_tm_win = np.empty(shape=0)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "calculated TM_WIN statistics",
        file=logfile,
    )


# fuse the summarizing stats list into a 2d-np.array, each row containing the
# summarizing stats of a single simulated tree sequence
dataframe_sumstats = np.concatenate(
    (dataframe_sumstats_sfs, dataframe_sumstats_ld, dataframe_sumstats_tm_win), axis=1
)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("calculated summarizing statisitcs", file=logfile)


# save treeseq list
np.save(snakemake.output.npy, dataframe_sumstats)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved sumstats to {snakemake.output.npy}", file=logfile)


# save the names of the sumstats
sumstat_type = (
    ["sfs"] * dataframe_sumstats_sfs.shape[1]
    + ["ld"] * dataframe_sumstats_ld.shape[1]
    + ["tm_win"] * dataframe_sumstats_tm_win.shape[1]
)
sumstat_type_id = [
    str(element)
    for element in list(range(dataframe_sumstats_sfs.shape[1]))
    + list(range(dataframe_sumstats_ld.shape[1]))
    + list(range(dataframe_sumstats_tm_win.shape[1]))
]

sumstat_names = np.array(list(map("_".join, zip(sumstat_type, sumstat_type_id))))


assert (
    len(sumstat_names) == dataframe_sumstats.shape[1]
), "name vector for summary statistics has the wrong dimension"


np.save(snakemake.output.sumstat_count, sumstat_names)
