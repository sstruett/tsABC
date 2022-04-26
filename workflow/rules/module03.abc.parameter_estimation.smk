"""
Here, we estimate the parameters for the performance analysis
"""


rule parameter_estimation:
    output:
        estims="results/abc/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.RDS",
    input:
        sumstats="results/abc/transformation/statcomp_{statcomposition}..simulations_sumstats.txt",
        podstats="results/abc/transformation/statcomp_{statcomposition}..pods_podstats.txt",
    log:
        log1="logs/module03/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
    script:
        "../scripts/02.parameter_estimation.R"


rule parameter_estimation_masked:
    output:
        estims="results/abc/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.masked.RDS",
    input:
        sumstats="results/abc/transformation/statcomp_{statcomposition}.simulations.sumstats.masked.txt",
        podstats="results/abc/transformation/statcomp_{statcomposition}.pods.podstats.masked.txt",
    log:
        log1="logs/module03/parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
    script:
        "../scripts/02.parameter_estimation.R"


rule aggregate_parameter_estimation:
    output:
        estims="results/abc/parameter_estimation.RDS",
    input:
        estims=expand(
            rules.parameter_estimation.output.estims,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module03/aggregate_parameter_estimation/parameter_estimation.log",
    run:
        sys.exit("#" * 600 + "inside aggregate_parameter_estimation\n" + "")


rule aggregate_parameter_estimation_masked:
    output:
        estims="results/abc/parameter_estimation_masked.RDS",
    input:
        estims=expand(
            rules.parameter_estimation_masked.output.estims,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module03/aggregate_parameter_estimation_masked/parameter_estimation_masked.log",
    run:
        sys.exit("#" * 600 + "inside aggregate_parameter_estimation_masked\n" + "")


rule plot_pod_sumstats:
    output:
        pdf="results/abc/pod_param_correlation.pdf",
    input:
        transformed=expand(
            "results/abc/transformation/statcomp_{statcomposition}..pods_podstats.txt",
            statcomposition=wildcards_statcomposition(config),
        ),
