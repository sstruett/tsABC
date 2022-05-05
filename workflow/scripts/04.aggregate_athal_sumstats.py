"""Aggregate podstats

Here, we aggregate the summarizing statistics of the pseudo-observed datasets
into a data frame. We provide the parameters as well as the pod ID as separate
columns do the data frame.
"""

import re
import datetime
import numpy as np
import pandas as pd


# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "start aggregating observed sumstats",
        file=logfile,
    )


# read dataframes and obtain popid from filenames
obsstats_list = []
identifier = []
for obsstats_filename in snakemake.input.observations:
    # obtain population identifier
    split = re.split(r"_|\.|\/", obsstats_filename)
    popid = split[np.where(np.array(split) == "population")[0].max() + 1]

    # read in data
    pddf = pd.read_feather(obsstats_filename)

    # create the map
    identifier.extend([popid] * len(pddf.index))

    obsstats_list.append(pddf)


pd.concat(obsstats_list, ignore_index=True).to_feather(
    snakemake.output.sumstats, compression="lz4"
)


with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved dataframe to {snakemake.output.sumstats}", file=logfile)


pd.DataFrame(identifier, columns=["population"]).to_feather(
    snakemake.output.identifier, compression="lz4"
)


with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved identifier to {snakemake.output.identifier}", file=logfile)
