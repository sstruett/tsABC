"""
Visualisation of the performance analysis and Arabidopsis thaliana

Although some exploratory and diagnostic plots were plottet during the pipeline,
this is the module to visualise the results into a somewhat more final state.
Here, we plot the paramter estimation, model choice and some investigatives
about parameters and pods.
"""


rule visualize_podstats:
    output:
        pdf="results/plots/abc/parameter.pdf",
    input:
        transformed=expand(
            "results/abc/transformation/statcomp_{statcomposition}..pods_podstats.txt",
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module05/visualize_podstats.log",
    conda:
        "config/env.yaml"
    script:
        "../scripts/05.visualize_podstats.R"


rule visualize_podstats_masked:
    output:
        pdf="results/plots/abc/parameter.masked.pdf",
    input:
        transformed=expand(
            "results/abc/transformation/statcomp_{statcomposition}..pods_podstats.masked.txt",
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module05/visualize_podstats_masked.log",
    conda:
        "config/env.yaml"
    script:
        "../scripts/05.visualize_podstats.R"


rule visualize_model_choice:
    output:
        pdf="results/plots/abc/model_choice.pdf",
    input:
        model_choice="results/abc/model_choice.RDS",
    log:
        log1="logs/module05/visualize_model_choice.log",
    conda:
        "config/env.yaml"
    params:
        # param no 4 (zero-based) must be t_sigma
        tsigma_per_podid=config["ABC"]["performance"]["pods"][3],
    script:
        "../scripts/05.visualize_model_choice.R"


rule visualize_model_choice_masked:
    output:
        pdf="results/plots/abc/model_choice.masked.pdf",
    input:
        model_choice="results/abc/model_choice_masked.RDS",
    log:
        log1="logs/module05/visualize_model_choice_masked.log",
    conda:
        "config/env.yaml"
    params:
        # param no 4 (zero-based) must be t_sigma
        tsigma_per_podid=config["ABC"]["performance"]["pods"][3],
    script:
        "../scripts/05.visualize_model_choice.R"


rule visualize_parameter_estimation:
    output:
        pdf="results/plots/abc/parameter_estimation.pdf",
    input:
        estimations="results/abc/parameter_estimation.RDS",
    log:
        log1="logs/module05/visualize_paremter_estimation.log",
    conda:
        "config/env.yaml"
    params:
        # param no 4 (zero-based) must be t_sigma
        tsigma_per_podid=config["ABC"]["performance"]["pods"][3],
        # which average interquantile ranges to calculate
        iqr=(0.99, 0.95, 0.9, 0.8, 0.5, 0.25, 0.1, 0),
    script:
        "../scripts/05.visualize_parameter_estimation.R"


rule visualize_parameter_estimation_masked:
    output:
        pdf="results/plots/abc/parameter_estimation.masked.pdf",
    input:
        estimations="results/abc/parameter_estimation_masked.RDS",
    log:
        log1="logs/module05/visualize_parameter_estimation_masked.log",
    conda:
        "config/env.yaml"
    params:
        # param no 4 (zero-based) must be t_sigma
        tsigma_per_podid=config["ABC"]["performance"]["pods"][3],
        # which average interquantile ranges to calculate
        iqr=(0.99, 0.95, 0.9, 0.8, 0.5, 0.25, 0.1, 0),
    script:
        "../scripts/05.visualize_parameter_estimation.R"
