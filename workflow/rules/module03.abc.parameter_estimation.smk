"""
Here, we estimate the parameters for the performance analysis
"""


rule parameter_estimation:
    output:
        estims="results/abc/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}...RDS",
        postplots=directory("results/abc/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}/"),
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
        statcompositions=config["ABC"]["sumstats"],
        statcomposition_no=lambda wildcards: int(float(wildcards.statcomposition)),
        regression="loclinear",  # r-abc package
    script:
        "../scripts/03.parameter_estimation.R"


rule parameter_estimation_masked:
    output:
        estims="results/abc/parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}...RDS",
        postplots=directory("results/abc/parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}/"),
    input:
        sumstats="results/abc/transformation/statcomp_{statcomposition}..simulations_sumstats.masked.txt",
        podstats="results/abc/transformation/statcomp_{statcomposition}..pods_podstats.masked.txt",
    log:
        log1="logs/module03/parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
        statcompositions=config["ABC"]["sumstats"],
        statcomposition_no=lambda wildcards: int(float(wildcards.statcomposition)),
        regression="loclinear",  # r-abc package
    script:
        "../scripts/03.parameter_estimation.R"


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
    script:
        "../scripts/03.aggregate_parameter_estimation.R"


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
    script:
        "../scripts/03.aggregate_parameter_estimation.R"


rule plot_pod_sumstats:
    output:
        pdf="results/abc/pod_param_correlation.pdf",
    input:
        transformed=expand(
            "results/abc/transformation/statcomp_{statcomposition}..pods_podstats.txt",
            statcomposition=wildcards_statcomposition(config),
        ),
