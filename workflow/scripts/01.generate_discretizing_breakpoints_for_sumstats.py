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
        # Loading list of tree sequences
        with open(snakemake.input.tsl[0], "rb") as tsl_file:  
            tsl = pickle.load(tsl_file)

        # random number generator to choose the haplotype per individual
        rng = np.random.default_rng(snakemake.params.seed)

        # calculate the data-based breakspoints
        breaks = pyfuncs.find_breakpoints_for_TM_WIN(
            tsl = tsl,
            specs = snakemake.config["ABC"]["sumstats_specs"]["TM_WIN"],
            rng = rng,
            log = snakemake.log.log1
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
