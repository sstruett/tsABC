# tsABC

## IMPORTANT
 + Why is there no LD stats in statcomp 0

## still to implement
 + calculation of masked sumstats
 + model choice (reduce dimensionality)
 + model choice for masked sumstats
 + parameter estimate (use reduced dimensionality provided in previous module)
 + parameter estimate for masked sumstats
 + update all R scripts in module 02; rfunctions R and the scripts that use them
 + rewrite the transformation part to use a single rule for the transformation

 + if no mask is provided, create a rule that makes a random mask, e.g. remove 30% of data on 50 kb stretches; e.g. put into readme how to generate this random mask; $ snakemake -j 1 results/mask/generate_random_mask.touch

## Possible improvement

 + for sims I average over loci in the aggregation, but for pods I do in the calculation of podstats script

# Updating files

## Workflow
 + snakefmt, $ for f in workflow/rules/*.smk; do echo $f; snakefmt $f; sleep 1; done; sleep 1; f="workflow/Snakefile"; echo $f; snakefmt $f
 + snakemake --lint

## Python scripts
 + black
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
