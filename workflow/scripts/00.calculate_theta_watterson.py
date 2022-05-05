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
import pickle
import pyfuncs  # from file


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        f"start calculating theta Watterson with mode: {snakemake.wildcards.mode}",
        file=logfile,
    )


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

        # log
        with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                f"chropped region {chromid} of "
                + f"""{len(snakemake.config["ABC"]["athaliana"]["observations"][
                    "treeseq_1001"]["chosen_region"])}""",
                file=logfile,
            )

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
                snakemake.config["ABC"]["athaliana"]["observations"]["num_observations"]
            )
        ),
        "nsam": int(float(snakemake.config["ABC"]["simulations"]["nsam"])),
    }
    tsl = pyfuncs.create_subsets_from_treeseqlist(
        treeseq_list, specs, rng, snakemake.log.log1
    )

elif snakemake.wildcards.mode == "genome":
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
    for chromid, start, stop in snakemake.config["ABC"]["athaliana"]["observations"][
        "treeseq_1001"
    ]["whole_gemome_approach"]:
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
        treeseq_list.append(
            (chromid, treeseq_athal_population.keep_intervals([(start, stop)]).trim())
        )

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
        print(
            "created region treeseq list (1 per region of each chromosome)",
            file=logfile,
        )

    # create subsample from treesequence
    specs = {
        "num_observations": int(
            float(
                snakemake.config["ABC"]["athaliana"]["observations"]["treeseq_1001"][
                    "whole_genome_approach_num_observations"
                ]
            )
        ),
        "nsam": int(float(snakemake.config["ABC"]["simulations"]["nsam"])),
    }

    # remove chromosome from treeseq_list
    treeseq_list = [treeseq for _, treeseq in treeseq_list]

    tsl = pyfuncs.create_subsets_from_treeseqlist(
        treeseq_list, specs, rng, snakemake.log.log1
    )

elif snakemake.wildcards.mode == "pod":
    # Loading list of tree sequences
    with open(snakemake.input.tsl_pod[0], "rb") as tsl_file:
        tsl = np.array(pickle.load(tsl_file), dtype=object)

else:
    assert False, "unknown mode to calclate theta Watterson"


# calculate theta_watterson per treeseq
theta_watterson = []
for treeseq in tsl.flatten():
    theta_watterson.append(
        treeseq.segregating_sites(span_normalise=True)
        / sum([1 / i for i in range(1, treeseq.num_samples)])
    )

theta_watterson = np.array(theta_watterson).mean()

# calculate expected 2N
two_N_zero = round(
    theta_watterson / (2 * float(snakemake.config["ABC"]["simulations"]["mutrate"]))
)


print("\n" + "_" * 80, file=sys.stderr)
print(f"theta_watterson\t{theta_watterson}", file=sys.stderr)
print(f"two_N_zero\t{two_N_zero}", file=sys.stderr)
print("=" * 80 + "\n", file=sys.stderr)

# print to output file
with open(snakemake.output.txt, "w", encoding="utf-8") as outfile:
    print(f"theta_watterson\t\t{theta_watterson}", file=outfile)
    print(f"two_N_zero\t\t{two_N_zero}", file=outfile)


# log the results
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"theta_watterson\t{theta_watterson}", file=logfile)
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"two_N_zero\t{two_N_zero}", file=logfile)
