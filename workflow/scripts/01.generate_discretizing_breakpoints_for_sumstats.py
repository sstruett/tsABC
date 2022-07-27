"""
Generate discretization boundaries for the summarizing statistics that need that

This is a rule that provides a suggestion. However, any other suggested breaks
may be provided instead of running this rule
"""


import datetime
import itertools
import numpy as np
import warnings
import pickle
import tskit
import json
import pyfuncs  # from file


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "start creating breaks for discretization of summarizing statistics",
        file=logfile,
    )


# get unique sumstats from sumstat sets
unique_sumstats = set(
    itertools.chain.from_iterable(
        [sumstat_set.split("/") for sumstat_set in snakemake.config["ABC"]["sumstats"]]
    )
)


# random number generator to choose the haplotype per individual and subsampling
rng = np.random.default_rng(snakemake.params.seed)


# loop through sumstats
for sumstat in unique_sumstats:
    breaks = None
    if sumstat == "SFS":
        breaks = np.empty(shape=0)  # we do not discretize SFS
    elif sumstat == "LD":
        breaks = np.array(
            [
                float(this_break)
                for this_break in snakemake.config["ABC"]["sumstats_specs"]["LD"][
                    "breaks"
                ]
            ]
        )
    elif sumstat == "TM_WIN":
        if snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]["breaks_mode"] == "pod":
            # Loading list of tree sequences
            with open(snakemake.input.tsl[0], "rb") as tsl_file:
                tsl = pickle.load(tsl_file)

            # calculate the data-based breakspoints
            breaks = pyfuncs.find_breakpoints_for_TM_WIN(
                tsl=tsl,
                specs=snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"],
                rng=rng,
                log=snakemake.log.log1,
            )
        elif (
            snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]["breaks_mode"]
            == "athal"
        ):
            # read tree sequence
            treeseq_athal = tskit.load(snakemake.input.tsl[0])

            # read sample names of first population
            sample_names = np.loadtxt(snakemake.input.tsl[1], dtype=str)

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

            # calculate the data-based breakspoints
            breaks = pyfuncs.find_breakpoints_for_TM_WIN(
                tsl=tsl,
                specs=snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"],
                rng=rng,
                log=snakemake.log.log1,
            )
        elif (
            snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]["breaks_mode"]
            == "expected"
        ):
            discretized_times = pyfuncs.discretized_times(
                n=snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]["classes"], M=2
            )
            mutrate = float(snakemake.config["ABC"]["simulations"]["mutrate"])
            two_N_zero = int(
                float(
                    snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"][
                        "two_N_zero_if_expected"
                    ]
                )
            )
            window_size = int(
                float(snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]["winsize"])
            )
            breaks = pyfuncs.snp_freq_from_times(
                discretized_times, two_N_zero, mutrate, window_size
            )
        else:
            sys.exit(
                "#" * 600
                + " inside generate_discretizing_breakpoints_for_sumstats\n"
                + f'your breaks_mode is maldefined: {snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"]["breaks_mode"]}\n'
                + "YOU SHOULD NEVER REACH HERE!"
            )

    # test if they belong to the sumstats that are implemented, may need to
    # co-check with the config check rule in module_00
    if breaks is None:
        warnings.warn(
            "".join(
                [
                    "some sumstats are not taken into account",
                    " for the creation of discretization breaks",
                ]
            )
        )

    outfile_name = f"resources/discretization/sumstat_{sumstat}.npytxt"
    np.savetxt(outfile_name, breaks)

    # log
    with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(f"created breaks for discretization of {sumstat}", file=logfile)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("created breaks for discretization of summarizing statistics", file=logfile)
