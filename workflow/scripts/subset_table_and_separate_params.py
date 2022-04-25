"""
reads transformer file, sumstat file and subsets the sumstats to the one chosen
for the transformer, add the params also to that file and write it down into a 
new unzipped .txt file
"""

import sys
import pandas as pd

transformer_file = sys.argv[1]
sumstat_file = sys.argv[2]
outfile_prefix = sys.argv[3]

# define outfile
outfile_sumstat = outfile_prefix + ".sumstat"

# read files
transformer = pd.read_csv(transformer_file, sep="\t", header=None)

if not sumstat_file.endswith(".txt"):
    sumstat = pd.read_feather(sumstat_file)
else:
    sumstat = pd.read_csv(sumstat_file, sep="\t")

# subset columns to wanted stats
colnames_to_be_transformed = [*transformer[0]]
sumstats_of_interest = sumstat[colnames_to_be_transformed]

# extract params
params_all = sumstat.filter(regex="^param_")

# finalize data frame, params + sumstats
subsetted_table = pd.concat([params_all, sumstats_of_interest], axis=1)

# write to file
subsetted_table.to_csv(outfile_sumstat, sep="\t", index=False)
