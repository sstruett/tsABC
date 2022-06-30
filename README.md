# tsABC

## How to run
 + provide correct parameters, e. g. mutation rate and recombination rate and nsim, then run everything

## IMPORTANT
 + confirm the downsampling for model choice works correctly
 + describe that `$ snakemake -j 1 theta_watterson_genome` will calculate theta watterson on the whole region while `$ snakemake -j 1 theta_watterson_region` does it only on the provided regions of athal, while `snakemake -j 1 theta_watterson_pod` will uses 1000 simulated loci for the values of the first pod
 + make sure to simulate from_ts with the correct model

## still to implement
 + do_visualisation to be implemented into the pipeline
 + implement proper ylimits for parameter estimation plots; make sure to save the podid into the pod tables
 + write the two functions to run the 6 parameter model for A. thaliana
 + for module 02 change script names to include module number
 + add parameter check for the athaliana priors
 + complete localrules in module04, maybe check also for the other modules
 + add checking mask files for length of locus
 + add checking that provided regions for athal are same as length of simulation
 + make plotting in target_files() function become true
 + summarize the model choice with the classification; produce a table also: pods(rows) --> percentage of each support class
 + configuration sanity; check that number of athal regions and loci are the same
 + create a summarizing rule for athal, that works on the whole genome, e.g. only non-pericentromeric, e.g. compare pericentromeric with non-pericentromeric
 + check config also breaks_mode for tm_win: should allow for three modes: ["pod", "athal", "expected"]
 + implement the whole regions into the whole pipeline and target rule
 + implement the plot of the sumstats of athal into the pipeline
 + implement the same way to plot the pars for the pods
 + make regression type being in the config file not only in the params
 + make all R scripts print to log file not stdout
 + separate module00.common.smk to module00.common.rules.smk and module00.common.func.smk
 + finsish visualisation
 + add athal visualisation
 + think about how to visua
 + provide different modes for LD breaks

 + if no mask is provided, create a rule that makes a random mask, e.g. remove 30% of data on 50 kb stretches; e.g. put into readme how to generate this random mask; $ snakemake -j 1 results/mask/generate_random_mask.touch
 + update all R scripts in module 02; rfunctions R and the scripts that use them
 + implement to observe the data on a much larger region of thaliana and make the estimate on the unmasked simulations
 + make calculate sumstats using the same script, no matter masked, nomask or proposed vs alternative model

## Suggested development improvement
 + make a pure piecewise-constant demography estimate; then use the theta over time to simulate the alternative model, while a pure change in rho according to the rescaling through selfing to make an estimate of the change in recombination rate, but use the change in selfing instead

# Updating files

## Workflow
 + snakefmt, e. g. $ for f in workflow/rules/*.smk; do echo $f; snakefmt $f; sleep 1; done; sleep 1; f="workflow/Snakefile"; echo $f; snakefmt $f
 + snakemake --lint

## Python scripts
 + black; e. g. $ for f in workflow/scripts/*.py; do echo $f; black $f ; done
 + pylint

## Yaml
 + ?


## Breaks that have been used for tsabc3
[0, 1.6413854204016185, 3.371536501050441, 5.200605743928798, 7.140593642054711, 9.205826318456989, 11.413598206039438, 13.785053314958535, 16.346419960511703, 19.130784024179853, 22.18070977791825, 25.552246278968695, 29.32130341997296, 33.59430798395769, 38.52712973842995, 44.3614195558365, 51.502013197891216, 60.707839516348194, 73.68272297580947, 95.86343275372768, inf]


                            {simorpod}     {dataset}

            "results/abc/   simulations/   sumstats                       .feather"
            "results/abc/   simulations/   sumstats.masked                .feather"
            "results/abc/   simulations/   alternative.sumstats           .feather"
            "results/abc/   simulations/   alternative.sumstats.masked    .feather"
            "results/abc/   pods/          podstats                       .feather"
            "results/abc/   pods/          podstats.masked                .feather"
