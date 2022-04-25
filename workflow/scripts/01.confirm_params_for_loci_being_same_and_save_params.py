"""Quick check that parameters of the different loci are the same

This should usually work, because the rng for drawing the parameters uses the
random seed of determined for the first locus. The random seeds for the
simulations are (pseudo-) independent.
"""


import datetime
import numpy as np


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "start confirming the parameters of the independent loci being equal",
        file=logfile,
    )


# read params
params_per_locus = [np.load(params_file) for params_file in snakemake.input.params]
LOCUS_INDEX = 0
for param_index in range(1, len(params_per_locus)):
    params1, params2 = params_per_locus[param_index - 1 : param_index + 1]
    assert (
        params1 == params2
    ).all(), "Parameters for the different loci are not the same!"
    with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(f"{LOCUS_INDEX} ?== {LOCUS_INDEX+1}: True", file=logfile)
        LOCUS_INDEX += 1


with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "confirmed the parameters of the independent loci being equal",
        file=logfile,
    )


# save one copy of the params into a .npy file
np.save(
    snakemake.output.params,
    # params1 is the pre-last np.array from the read files, but as all of them
    # are checked for equality, it does not matter which one to save
    params1,
)


with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved params to .npy file: {snakemake.output.params}", file=logfile)
