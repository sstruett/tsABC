"""
simulate the treesequences under the six parameter model: A single change in
selfing rate and a single change in population size
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
    seed = np.load(snakemake.input.npy)[simid][int(float(snakemake.wildcards.locid))]
    seed_params = np.load(snakemake.input.npy)[simid][
        0
    ]  # the parameters for the independent loci must be the same
    rng = np.random.default_rng(seed)
    rng_params = np.random.default_rng(seed_params)

    # draw the parameters of the model
    params = [
        pyfuncs.draw_parameter_from_prior(prior_definition, rng_params)
        for prior_definition in snakemake.params.model["priors"]
    ]
    del (
        seed_params,
        rng_params,
    )  # only the parameter drawing relies on the same seed, the simulations must be independent

    # simulate
    ts = pyfuncs.simulate_treesequence_under_model(
        params, snakemake.params, rng, snakemake.log.log1
    )

    return ts, params


# define the loop for the clustering of consequential simulations
start = int(float(snakemake.wildcards.simid))
end = min(
    start + int(float(snakemake.config["ABC"]["jobluster_simulations"])),
    int(float(snakemake.config["ABC"]["simulations"]["nsim"])),
)

assert start < end, "simid (Simulation IDs) maldefined, you will not simulate anything"

# run the simulations
treesequence_list = []
parameters_list = []
for simid in range(start, end):
    # log
    with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(f"creating treesequence no {simid}", file=logfile)

    treeseq, params = main(simid)
    treesequence_list.append(treeseq)
    parameters_list.append(params)


# save treeseq list
with open(snakemake.output.tsl, "wb") as tsl_file:
    pickle.dump(treesequence_list, tsl_file)


# save parameters list
np.save(snakemake.output.params, np.array(parameters_list))


# debug
sys.exit("#"*600 + " what is happening here? simulate six par model")

