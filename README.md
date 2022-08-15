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


## still to implement
 + describe that `$ snakemake -j 1 theta_watterson_genome` will calculate theta watterson on the whole region while `$ snakemake -j 1 theta_watterson_region` does it only on the provided regions of athal, while `snakemake -j 1 theta_watterson_pod` will uses 1000 simulated loci for the values of the first pod
 + implement proper ylimits for parameter estimation plots; make sure to save the podid into the pod tables
 + for module 02 change script names to include module number
 + add parameter check for the athaliana priors
 + add checking mask files for length of locus
 + make plotting in target_files() function become true
 + configuration sanity; check that number of athal regions and loci are the same
 + make regression type being in the config file not only in the params
 + make all R scripts print to log file not stdout
 + provide different modes for LD breaks

 + probably provide random masking if no masking file is provided
 + update all R scripts in module 02; rfunctions R and the scripts that use them

## Suggested development improvement
 + consequential ABC with increasing complexity of models
