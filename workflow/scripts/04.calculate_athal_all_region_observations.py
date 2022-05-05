"""
caluclate summarizing statistics on the inferred treesequences of the 1001
genomes. This script does following steps:
  + loading and subsetting treeseq to the population;
  + extracting one treesequence per chosen region for each chromosome;
  + masking if asked for;
  + create a maximum number of subsets of the configured sample size;
  + calculate the summary statistics on each of these subsets;
  + average the summary statistics over the different chromosomal regions
"""


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


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "start calculating observed summarizing statistics on Arabidopsis thaliana data",
        file=logfile,
    )


# read sample names
sample_names = np.loadtxt(snakemake.input.sample_list, dtype=str)

# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        f"read sample names from {snakemake.input.sample_list}",
        file=logfile,
    )


# read tree sequence
treeseq_athal = tskit.load(snakemake.input.athal_treeseq)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        f"read treeseq of Arabidopsis thaliana from '{snakemake.input.athal_treeseq}'",
        file=logfile,
    )


# find the node ids for the sample of the population
population_sample = []
for individual in treeseq_athal.individuals():
    if str(json.loads(individual.metadata)["id"]) in sample_names:
        population_sample.extend(individual.nodes)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("created the node id sample set for the given population", file=logfile)


# sample treeseq to provided samples
treeseq_athal_population = treeseq_athal.simplify(samples=population_sample)
del treeseq_athal


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("simplified the tree sequence to the given population", file=logfile)


# get the chromosomal regions from the config file
chromosome_regions = []
for chromid, start, stop in snakemake.config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["whole_gemome_approach"]:
    start += snakemake.params.chrom_multiplier * chromid
    stop += snakemake.params.chrom_multiplier * chromid
    chromosome_regions.append((start, stop))


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("prepared regions", file=logfile)


# chop down to regions
treeseq_list = []
for start, stop in chromosome_regions:
    # calculate chromid from position in original treeseq
    chromid = int(start / snakemake.params.chrom_multiplier)

    # add tuple with chromosome/treeseq, as chromid is needed for masking
    treeseq_list.append((chromid, treeseq_athal_population.keep_intervals([(start, stop)]).trim()))

    # log
    with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(
            f"prepared regions for chromsome {chromid + 1} of {len(chromosome_regions)}",
            file=logfile,
        )

del treeseq_athal_population


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("created region treeseq list (1 per region of each chromosome)", file=logfile)


# mask if asked; make sure the parameter is a boolean and then mask or not
do_mask = pyfuncs.check_masking_parameter(snakemake.params.masked)

# read and prepare mask files
if do_mask:
    # provide mask in parallel to the regions of the treeseq_list
    mask = []
    for treeseqid, (chromid, _) in enumerate(treeseq_list):
        this_mask = np.loadtxt(snakemake.input.maskfiles[chromid - 1], dtype=int)
        region_start, region_end = chromosome_regions[treeseqid]
        region_start -= chromid * snakemake.params.chrom_multiplier
        region_end -= chromid * snakemake.params.chrom_multiplier
        region_mask = this_mask[
            (region_start <= this_mask[:, 0]) & (region_end > this_mask[:, 1])
        ]
        this_mask = region_mask - region_start
        del region_mask

        # filter mask for disjoint intervals
        mask.append(
            pyfuncs.filter_mask_for_disjoint_intervals(
                this_mask, log=snakemake.log.log1
            )
        )

    # remove chromid from treeseq_list
    treeseq_list = [treeseq for _, treeseq in treeseq_list]


    # log
    with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(
            f"prepared coordinates to mask the exons for {len(mask)} regions",
            file=logfile,
        )

    treeseq_list_masked = []
    for treeid, (treeseq_chrom, mask_chrom) in enumerate(
        zip(treeseq_list, mask), start=1
    ):
        treeseq_list_masked.append(
            treeseq_chrom.delete_intervals(
                mask_chrom, simplify=True, record_provenance=True
            )
        )

        # log
        with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                f"masked {treeid} of {len(mask)} chromosomes",
                file=logfile,
            )

    treeseq_list = treeseq_list_masked
    del treeseq_list_masked
else:
    # remove chromid from treeseq_list
    treeseq_list = [treeseq for _, treeseq in treeseq_list]


    # log, create log file
    with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("will not mask for exons", file=logfile)


# create subsample from treesequence
rng = np.random.default_rng(snakemake.params.seed)
specs = {
    "num_observations": int(
        float(snakemake.config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["whole_genome_approach_num_observations"])
    ),
    "nsam": int(float(snakemake.config["ABC"]["simulations"]["nsam"])),
}
tsl = pyfuncs.create_subsets_from_treeseqlist(
    treeseq_list, specs, rng, snakemake.log.log1
)


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

    dataframe_sumstats_sfs = np.empty(tsl.shape, dtype=tskit.TreeSequence)
    for my_id, (treeid, treeseq) in enumerate(np.ndenumerate(tsl), start=1):
        dataframe_sumstats_sfs[treeid] = treeseq.allele_frequency_spectrum(
            sample_sets=None,
            windows=None,
            mode="site",
            span_normalise=False,
            polarised=snakemake.config["ABC"]["sumstats_specs"]["SFS"]["polarised"],
        )

        # log
        if not my_id % 100:
            with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"calculating SFS: {my_id} of {tsl.size}",
                    file=logfile,
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
    dataframe_sumstats_sfs = np.empty(shape=(0, 0))


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
    for my_id, (treeid, treeseq) in enumerate(np.ndenumerate(tsl), start=1):
        dataframe_sumstats_ld[treeid] = pyfuncs.calculate_ld(
            treeseq, specs, breaks, rng, snakemake.log.log1
        )

        # log
        if not my_id % 5:
            with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"calculating LD: {my_id} of {tsl.size}",
                    file=logfile,
                )

    # get average of ld over the different loci
    dataframe_sumstats_ld = dataframe_sumstats_ld.mean(axis=1)

    dataframe_sumstats_ld = np.concatenate(dataframe_sumstats_ld, axis=0).reshape(
        (len(dataframe_sumstats_ld), dataframe_sumstats_ld[0].shape[0])
    )
else:
    dataframe_sumstats_ld = np.empty(shape=(0, 0))


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
    for my_id, (treeid, treeseq) in enumerate(np.ndenumerate(tsl), start=1):
        dataframe_sumstats_tm_win[treeid] = pyfuncs.calculate_tm_win(
            treeseq, specs, breaks_tm_win
        )

        # log
        if not my_id % 100:
            with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"calculating TM_WIN: {my_id} of {tsl.size}",
                    file=logfile,
                )

    # get average of tm_win over the different loci
    dataframe_sumstats_tm_win = dataframe_sumstats_tm_win.mean(axis=1)
    dataframe_sumstats_tm_win = np.concatenate(
        dataframe_sumstats_tm_win, axis=0
    ).reshape((len(dataframe_sumstats_tm_win), dataframe_sumstats_tm_win[0].shape[0]))
else:
    dataframe_sumstats_tm_win = np.empty(shape=(0, 0))


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
    [
        npy
        for npy in (
            dataframe_sumstats_sfs,
            dataframe_sumstats_ld,
            dataframe_sumstats_tm_win,
        )
        if npy.size != 0
    ],
    axis=1,
)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("calculated summarizing statisitcs", file=logfile)


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


pd.DataFrame(dataframe_sumstats, columns=sumstat_names).to_feather(
    snakemake.output.sumstats, compression="lz4"
)

# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved sumstats to {snakemake.output.sumstats}", file=logfile)
