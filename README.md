# tsABC

## How to run
 + provide correct parameters, e. g. mutation rate and recombination rate and nsim, then run everything

## IMPORTANT
 + Why is there no LD stats in statcomp 0, maybe use more SNPs to calculate all classes, but make sure that NaN is not a problem
 + confirm the downsampling for model choice works correctly

## still to implement
 + parameter estimate
 + parameter estimate for masked sumstats (should work using the same script)
 + write the two functions to run the 6 parameter model for A. thaliana
 + Obtain stats from A. thaliana
 + for module 02 change script names to include module number
 + add parameter check for the athaliana priors
 + complete localrules in module04, maybe check also for the other modules
 + add checking mask files for length of locus
 + add checking that provided regions for athal are same as length of simulation
 + make plotting in target_files() function become true
 + summarize the model choice with the classification; produce a table also: pods(rows) --> percentage of each support class
 + configuration sanity; check that number of athal regions and loci are the same
 + create a summarizing rule for athal, that works on the whole genome, e.g. only non-pericentromeric, e.g. compare pericentromeric with non-pericentromeric

 + if no mask is provided, create a rule that makes a random mask, e.g. remove 30% of data on 50 kb stretches; e.g. put into readme how to generate this random mask; $ snakemake -j 1 results/mask/generate_random_mask.touch
 + update all R scripts in module 02; rfunctions R and the scripts that use them
 + implement to observe the data on a much larger region of thaliana and make the estimate on the unmasked simulations

## Possible improvement

 + For the sims I average over loci in the aggregation, but for pods I do in the calculation of podstats script, more consistency needed

# Updating files

## Workflow
 + snakefmt, e. g. $ for f in workflow/rules/*.smk; do echo $f; snakefmt $f; sleep 1; done; sleep 1; f="workflow/Snakefile"; echo $f; snakefmt $f
 + snakemake --lint

## Python scripts
 + black; e. g. $ for f in workflow/scripts/*.py; do echo $f; black $f ; done
 + pylint

## Yaml
 + ?



                            {simorpod}     {dataset}

            "results/abc/   simulations/   sumstats                       .feather"
            "results/abc/   simulations/   sumstats.masked                .feather"
            "results/abc/   simulations/   alternative.sumstats           .feather"
            "results/abc/   simulations/   alternative.sumstats.masked    .feather"
            "results/abc/   pods/          podstats                       .feather"
            "results/abc/   pods/          podstats.masked                .feather"
