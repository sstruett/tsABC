"""
simulate the treesequences
"""


import sys
import datetime
import pickle
import numpy as np
import tskit
import msprime
import pyfuncs  # from file


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("start simulating tree sequences", file=logfile)


def main(simid):
    """Heart of this script

    The main function will create a single tree sequence. The main function has
    to be repeatedly executed until the list of treesequences has been created.
    """
    # provide the correct seed to create the random number generator
    seed = np.load(snakemake.input.npy)[
        simid[0],
        simid[1],
        int(float(snakemake.wildcards.podid)),
    ]
    rng = np.random.default_rng(seed)

    # draw the parameters of the model
    params = (
        int(
            float(snakemake.params["params"][0][int(float(snakemake.wildcards.podid))])
        ),
        float(snakemake.params["params"][1][int(float(snakemake.wildcards.podid))]),
        float(snakemake.params["params"][2][int(float(snakemake.wildcards.podid))]),
        int(
            float(snakemake.params["params"][3][int(float(snakemake.wildcards.podid))])
        ),
    )

    print(params)
    sys.exit(
        "#" * 600
        + " inside 01.simulate_treeseq_pods.py, find the issue with the params"
    )

    # simulate
    ts = pyfuncs.simulate_treesequence_under_model(
        params, snakemake.params, rng, snakemake.log.log1
    )

    return ts


# run the simulations
pod_list_of_treesequence_list_per_locus = []
for repid in range(snakemake.params.nsim):
    treesequence_list = []
    for locid in range(snakemake.params.nloci):
        # log
        with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                f"creating treeseq for pod {snakemake.wildcards.podid} locus {locid} of rep {repid}",
                file=logfile,
            )

        treeseq = main((repid, locid))
        treesequence_list.append(treeseq)
    pod_list_of_treesequence_list_per_locus.append(treesequence_list)


# save treeseq list
with open(snakemake.output.tsl, "wb") as tsl_file:
    pickle.dump(pod_list_of_treesequence_list_per_locus, tsl_file)
