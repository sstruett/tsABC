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
        "start aggregating params and sumstats",
        file=logfile,
    )


# obtain simid and locid from filenames
npodid = len(snakemake.params.podid_wc)

# create data list to read and store by locus
read_data_list = [[] for _ in range(npodid)]
for npy in snakemake.input.npys:
    split = re.split(r"_|\.|\/", npy)
    podid = int(split[np.where(np.array(split) == "podid")[0].max() + 1])
    read_data_list[
        np.where(podid == np.array(snakemake.params.podid_wc))[0][0]
    ] = np.load(npy)
read_data = np.array(read_data_list, dtype=float)

with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("worked out sumstats", file=logfile)


# read sumstat names
sumstat_names_list = []
for name_file in snakemake.input.sumstat_names:
    sumstat_names_list.append(np.load(name_file))
    if len(sumstat_names_list) >= 2:
        assert all(
            sumstat_names_list[-2] == sumstat_names_list[-1]
        ), "sumstats are not congruent"
sumstat_names = sumstat_names_list[0]


# read params from config
npodid = len(snakemake.config["ABC"]["performance"]["pods"][0])
param_values = [[] for _ in range(npodid)]
for podid in range(npodid):
    for param in snakemake.config["ABC"]["performance"]["pods"]:
        param_values[podid].append(param[podid])
param_values = np.array(param_values).astype(float)


# concatenate and parse data into pandas dataframe; add the parameters
result_table_list = []
for podid in range(read_data.shape[0]):
    result_table_list.append(pd.DataFrame(data=read_data[podid], columns=sumstat_names))
    result_table_list[-1]["podid"] = podid

    parameter_columns_names = []
    for parid, param in enumerate(param_values[podid]):
        result_table_list[-1][f"param_{parid}"] = param
        parameter_columns_names.append(f"param_{parid}")

    # move parameter columns to front
    result_table_list[-1] = result_table_list[-1][
        parameter_columns_names
        + [
            col
            for col in result_table_list[-1].columns
            if col not in parameter_columns_names
        ]
    ]


pd.concat(result_table_list, ignore_index=True).to_feather(
    snakemake.output.sumstats, compression="lz4"
)

with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved dataframe to {snakemake.output.sumstats}", file=logfile)
