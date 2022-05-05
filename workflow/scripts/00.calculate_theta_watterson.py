"""Calculate an estimate of theta Watterson and two_N_zero

This is meant to be a proper estimate to use for the config file for proper
discretization of tm_win
"""


import sys
import datetime
import numpy as np
import msprime
import tskit
import json
import pyfuncs  # from file


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"start calculating theta Watterson with mode: {snakemake.wildcards.mode}", file=logfile)


# random seed generator subsetting and simulations
rng = np.random.default_rng(snakemake.params.seed)


if snakemake.wildcards.mode == "region":
    # read tree sequence
    treeseq_athal = tskit.load(snakemake.input.treeseq_athal)


    # read sample names of first population
    sample_names = np.loadtxt(snakemake.input.samples, dtype=str)


    # find the node ids for the sample of the population
    population_sample = []
    for individual in treeseq_athal.individuals():
        if str(json.loads(individual.metadata)["id"]) in sample_names:
            population_sample.extend(individual.nodes)

    # sample treeseq to provided samples
    treeseq_athal_population = treeseq_athal.simplify(samples=population_sample)
    del treeseq_athal

    # get the chromosomal regions from the config file
    chromosome_regions = []
    for chromid, (start, stop) in enumerate(
        snakemake.config["ABC"]["athaliana"]["observations"]["treeseq_1001"][
            "chosen_region"
        ],
        start=1,
    ):
        start += snakemake.params.chrom_multiplier * chromid
        stop += snakemake.params.chrom_multiplier * chromid
        chromosome_regions.append((start, stop))

        # log, create log file
        with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(f'chropped region {chromid} of {len(snakemake.config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["chosen_region"])}', file=logfile)

    # chop down to regions
    treeseq_list = []
    for chromid, (start, stop) in enumerate(chromosome_regions):
        treeseq_list.append(
            treeseq_athal_population.keep_intervals([(start, stop)]).trim()
        )
    del treeseq_athal_population

    # create subsample from treesequence
    specs = {
        "num_observations": int(
            float(
                snakemake.config["ABC"]["athaliana"]["observations"][
                    "num_observations"
                ]
            )
        ),
        "nsam": int(float(snakemake.config["ABC"]["simulations"]["nsam"])),
    }
    tsl = pyfuncs.create_subsets_from_treeseqlist(
        treeseq_list, specs, rng, snakemake.log.log1
    )

    print(tsl)
    print(tsl.shape)

    # calculate theta_watterson per treeseq
    theta_watterson = []
    for treeseq in tsl:
        theta_watterson.append(treeseq.segregating_sites()/sum([1/i for i in range(1, treeseq.num_samples)]))

    theta_watterson = np.array(theta_watterson).mean()
elif snakemake.wildcards.mode == "genome":
    sys.exit("implement mode 'genome'")
elif snakemake.wildcards.mode == "pod":
    sys.exit("implement mode 'pod'")
else:
    assert False, "unknown mode to calclate theta Watterson"

print(theta_watterson)