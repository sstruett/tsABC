"""
creating random seeds for simulations
"""


import datetime
import numpy as np


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("start creating random seeds for simulations", file=logfile)


# create new random seed generator
rng = np.random.default_rng(snakemake.params.seed)

# draw seeds and save as .npy to file
np.save(
    snakemake.output.npy,
    rng.integers(low=0, high=np.iinfo(int).max, size=snakemake.params.nseed),
)


# log
with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("created random seeds for simulations", file=logfile)
