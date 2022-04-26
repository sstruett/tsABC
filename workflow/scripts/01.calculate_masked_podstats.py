"""
caluclate summarizing statistics on the simulated treesequences
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
    print("start calculating summarizing statistics", file=logfile)


# Loading list of tree sequences
with open(snakemake.input.tsl, "rb") as tsl_file:
    tsl = np.array(pickle.load(tsl_file))


# sample only one haplotype per individual
rng = np.random.default_rng(snakemake.params.seed)
tsl_haploid = np.empty(tsl.shape, dtype=tskit.TreeSequence)
for treeid, treeseq in np.ndenumerate(tsl):
    sample_set = [
        this_individual.nodes[rng.integers(low=0, high=2)]
        for this_individual in treeseq.individuals()
    ]
    tsl_haploid[treeid] = treeseq.simplify(samples=sample_set)
del tsl
tsl = tsl_haploid
del tsl_haploid


# mask for regions provided
for treeseq in tsl:
    print(snakemake.input.mask)


    sys.exit("#" * 600 + "must implement the masking")


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("read and masked the tree sequences", file=logfile)


# read the breakpoint files
for breaks_filename in snakemake.input.breakpoints:
    if "SFS" in breaks_filename:
        with warnings.catch_warnings():  # will prevent warning if file is empty
            warnings.simplefilter("ignore")
            breaks_sfs = np.loadtxt(breaks_filename)
    elif "LD" in breaks_filename:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            breaks_ld = np.loadtxt(breaks_filename)
    elif "TM_WIN" in breaks_filename:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            breaks_tm_win = np.loadtxt(breaks_filename)

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

    dataframe_sumstats_sfs = np.empty(tsl.shape, dtype=tskit.TreeSequence)
    for treeid, treeseq in np.ndenumerate(tsl):
        dataframe_sumstats_sfs[treeid] = treeseq.allele_frequency_spectrum(
            sample_sets=None,
            windows=None,
            mode="site",
            span_normalise=False,
            polarised=snakemake.config["ABC"]["sumstats_specs"]["SFS"]["polarised"],
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

    # get average of sfs over the different loci
    dataframe_sumstats_sfs = dataframe_sumstats_sfs.mean(axis=1)

    dataframe_sumstats_sfs = np.concatenate(dataframe_sumstats_sfs, axis=0).reshape(
    (len(dataframe_sumstats_sfs), dataframe_sumstats_sfs[0].shape[0])
    )
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

    dataframe_sumstats_ld = np.empty(tsl.shape, dtype=tskit.TreeSequence)
    for treeid, treeseq in np.ndenumerate(tsl):
        dataframe_sumstats_ld[treeid] = pyfuncs.calculate_ld(
            treeseq, specs, breaks, rng, snakemake.log.log1
        )

    # get average of ld over the different loci
    dataframe_sumstats_ld = dataframe_sumstats_ld.mean(axis=1)

    dataframe_sumstats_ld = np.concatenate(dataframe_sumstats_ld, axis=0).reshape(
    (len(dataframe_sumstats_ld), dataframe_sumstats_ld[0].shape[0])
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

    dataframe_sumstats_tm_win = np.empty(tsl.shape, dtype=tskit.TreeSequence)
    for treeid, treeseq in np.ndenumerate(tsl):
        dataframe_sumstats_tm_win[treeid] = pyfuncs.calculate_tm_win(
            treeseq, specs, breaks_tm_win
        )

    # get average of tm_win over the different loci
    dataframe_sumstats_tm_win = dataframe_sumstats_tm_win.mean(axis=1)
    dataframe_sumstats_tm_win = np.concatenate(dataframe_sumstats_tm_win, axis=0).reshape(
    (len(dataframe_sumstats_tm_win), dataframe_sumstats_tm_win[0].shape[0])
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
