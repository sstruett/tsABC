import numpy as np
import pandas as pd
import re
import datetime

# log, create log file
with open(snakemake.log.log1, "w", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(
        "start aggregating params and sumstats",
        file=logfile,
    )


# obtain simid and locid from filenames
nloci = len(snakemake.params.locid_wc)
nsimid = len(snakemake.params.simid_wc)

# create data list to read and store by locus
read_data_list = [[[] for _ in range(nloci)] for _ in range(nsimid)]
for npy in snakemake.input.npys:
    split = re.split(r"_|\.|\/", npy)
    simid = int(split[np.where(np.array(split) == "sim")[0].max() + 1])
    locid = int(split[np.where(np.array(split) == "locus")[0].max() + 1])
    read_data_list[np.where(simid == np.array(snakemake.params.simid_wc))[0][0]][
        locid
    ] = np.load(npy)
read_data = np.array(read_data_list)


# take mean value over the independent loci
read_data = read_data.mean(axis=1)


# reshape the data into a 2d-array, because each input file provides a cluster
# batch of independent simulations
read_data = read_data.reshape(
    (read_data.shape[0] * read_data.shape[1], read_data.shape[2])
)

with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("worked out sumstats", file=logfile)

# read parameters
read_param_list = [[] for _ in range(nsimid)]
for npy in snakemake.input.params:
    split = re.split(r"_|\.|\/", npy)
    simid = int(split[np.where(np.array(split) == "sim")[0].max() + 1])
    read_param_list[
        np.where(simid == np.array(snakemake.params.simid_wc))[0][0]
    ] = np.load(npy)
read_param = np.array(read_param_list)
num_provided_params = read_param.shape[2]
read_param = read_param.reshape(
    (read_param.shape[0] * read_param.shape[1], read_param.shape[2])
)

with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("worked out params", file=logfile)

# concatenate parameters and sumstatas
result_table = np.concatenate((read_param, read_data), axis=1)

with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print("concatenated out params and sumstats", file=logfile)


# read sumstat names
sumstat_names_list = []
for name_file in snakemake.input.sumstat_names:
    sumstat_names_list.append(np.load(name_file))
    if len(sumstat_names_list) >= 2:
        assert all(
            sumstat_names_list[-2] == sumstat_names_list[-1]
        ), "sumstats are not congruent"
sumstat_names = sumstat_names_list[0]

# add param names
column_names = np.concatenate(
    (
        np.array([f"param_{parid}" for parid in range(num_provided_params)]),
        sumstat_names,
    ),
    axis=0,
)

result_table = pd.DataFrame(data=result_table, columns=column_names)
result_table.to_feather(snakemake.output.sumstats, compression="lz4")

with open(snakemake.log.log1, "a", encoding="utf-8") as logfile:
    print(datetime.datetime.now(), end="\t", file=logfile)
    print(f"saved dataframe to {snakemake.output.sumstats}", file=logfile)
