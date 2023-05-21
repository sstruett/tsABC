# tsABC


This git provides a snakemake workflow to run tsABC, the method introduced in our paper doi: https://doi.org/10.1101/2022.07.29.502030


The workflow requires some knowledge about snakemake (https://snakemake.readthedocs.io/en/stable/).

Theoretically, the whole pipeline could be run with $ snakemake -j 1
However, we suggest a more careful use with exact understanding of each step.

tsABC was first developed to run on A. thaliana, thus, the data application module is named "athal" (module04).


## How to run
 + provide correct parameters, e. g. mutation rate and recombination rate and nsim, then run everything

## IMPORTANT
 + The downsampling for model choice does NOT work correctly; there are two possible ways to handle it: First, one could simulate the same number of simulations for both proposed models. Second, the downsampling should be taken care manually in the according (or before) R script.


