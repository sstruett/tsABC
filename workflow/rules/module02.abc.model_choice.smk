"""
Here, we perform a model choice between the proposed and the confounding model
"""


import sys


rule find_pls:
    output:
        transformer="results/abc/transformation/transformer_{statcomposition}/Routput_{statcomposition}",
        rmse_plot="results/abc/transformation/transformer_{statcomposition}/RMSE_{statcomposition}.pdf",
    input:
        sumstats=rules.aggregate_sumstats.output.sumstats,
        script="workflow/scripts/find_pls_by_wanted_sumstats.R",
    log:
        log1="logs/module02/find_pls/transformed_{statcomposition}.log",
    conda:
        "config/env.yaml"
    shadow:
        "shallow"
    threads: 1
    params:
        num_pls_max=config["ABC"]["performance"]["pls_max"],
        sumstats_to_use=lambda wildcards: config["ABC"]["sumstats"][
            int(float(wildcards.statcomposition))
        ],
    shell:
        r"""
        # please leave threads 1, otherwise writing between cores seems to slow down the process
        #export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}

        which Rscript
        which R
        Rscript --vanilla {input.script} \
            {output.transformer} \
            {output.rmse_plot} \
            {input.sumstats} \
            {params.num_pls_max} \
            {wildcards.statcomposition} \
        "{params.sumstats_to_use}"

        """


rule find_pls_masked:
    output:
        transformer="results/abc/transformation/masked_transformer_{statcomposition}/Routput_{statcomposition}",
        rmse_plot="results/abc/transformation/masked_transformer_{statcomposition}/RMSE_{statcomposition}.pdf",
    input:
        sumstats=rules.aggregate_sumstats_masked.output.sumstats,
        script="workflow/scripts/find_pls_by_wanted_sumstats.R",
    log:
        log1="logs/module02/find_pls/transformed_{statcomposition}.log",
    conda:
        "config/env.yaml"
    shadow:
        "shallow"
    threads: 1
    params:
        num_pls_max=config["ABC"]["performance"]["pls_max"],
        sumstats_to_use=lambda wildcards: config["ABC"]["sumstats"][
            int(float(wildcards.statcomposition))
        ],
    shell:
        r"""
        # please leave threads 1, otherwise writing between cores seems to slow down the process
        #export OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads}

        which Rscript
        which R
        Rscript --vanilla {input.script} \
            {output.transformer} \
            {output.rmse_plot} \
            {input.sumstats} \
            {params.num_pls_max} \
            {wildcards.statcomposition} \
        "{params.sumstats_to_use}"

        """


rule transform_sumstats:
    output:
        sumstats="results/abc/transformation/statcomp_{statcomposition}..{simorpod}_{dataset}.txt",
    input:
        sumstats="results/abc/{simorpod}/{dataset}.feather",
        transformer=which_transformer,
        script="workflow/scripts/transformer",
        subsetter="workflow/scripts/subset_table_and_separate_params.py",
    log:
        log1="logs/module02/transform_sumstats/transformation_{statcomposition}.{simorpod}.{dataset}.log",
    conda:
        "config/env.yaml"
    shadow:
        "shallow"
    group:
        "transformation"
    threads: 1
    shell:
        r"""
        # This script transforms the summary statistics and the statistics of
        # the pseudo-observed datasets using the previously calculated pls; we
        # use the tools from the ABCtoolbox to do so.
        # The data has to be provided as .txt file for the transformer; this is
        # a pre-step to the actual calculation.


        # step 1: subset the stats/podstats and provide them as .txt
        python {input.subsetter} {input.transformer} {input.sumstats} subset.sumstat.table

        # log
        date > {log.log1}
        echo "subsetted and wrote sumstats and podstats as .txt" >> {log.log1}


        # step 2: transform the stats/podstats
        ./{input.script} {input.transformer} subset.sumstat.table.sumstat {output.sumstats} boxcox

        # log
        date >> {log.log1}
        echo "transformed sumstats and podstats" >> {log.log1}


        # give head of both files
        echo "\n\nsumstats' head" >> {log.log1}
        head {output.sumstats} >> {log.log1}

        """


rule model_choice:
    output:
        bayes_factors="results/abc/model_choice/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}...RDS",
        bayes_plot="results/abc/model_choice/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.pdf"
    input:
        sumstats="results/abc/transformation/statcomp_{statcomposition}..simulations_sumstats.txt",
        alternative_sumstats="results/abc/transformation/statcomp_{statcomposition}..simulations_alternative.sumstats.txt",
        podstats="results/abc/transformation/statcomp_{statcomposition}..pods_podstats.txt",
    log:
        log1="logs/module02/model_choice/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
    script:
        "../scripts/02.model_choice.R"


rule model_choice_masked:
    output:
        bayes_factors="results/abc/model_choice/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}...masked.RDS",
        bayes_plot="results/abc/model_choice/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.masked.pdf"
    input:
        sumstats="results/abc/transformation/statcomp_{statcomposition}..simulations_sumstats.masked.txt",
        alternative_sumstats="results/abc/transformation/statcomp_{statcomposition}..simulations_alternative.sumstats.masked.txt",
        podstats="results/abc/transformation/statcomp_{statcomposition}..pods_podstats.masked.txt",
    log:
        log1="logs/module02/model_choice_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
    script:
        "../scripts/02.model_choice.R"


rule aggregate_model_choice:
    output:
        bayes_factors="results/abc/model_choice.RDS",
    input:
        bayes_factors=expand(
            rules.model_choice.output.bayes_factors,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module02/aggregate_model_choice/model_choice.log",
    run:
        sys.exit("#" * 600 + "inside aggregate_model_choice\n" + "")


rule aggregate_model_choice_masked:
    output:
        bayes_factors="results/abc/model_choice_masked.RDS",
    input:
        bayes_factors=expand(
            rules.model_choice_masked.output.bayes_factors,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module02/aggregate_model_choice_masked/model_choice.log",
    run:
        sys.exit("#" * 600 + "inside aggregate_model_choice_masked\n" + "")
