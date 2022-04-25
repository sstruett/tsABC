"""
Here, we simulate the data under the designed model and calculate the
summarizing statistics. The final output of this module is a huge dataframe
to provide parameters and according summary statistics.
"""


import itertools


localrules:
    generate_discretizing_breakpoints_for_sumstats,
    confirm_params_for_loci_being_same_and_save_params,


rule create_abc_simulation_randints:
    output:
        npy=protected("resources/rand.int.npy"),
    log:
        log1="logs/module01/create_abc_simulation_randints.log",
    conda:
        "config/env.yaml"
    params:
        seed=int(float(config["ABC"]["seed"])),
        nseed=(
            int(float(config["ABC"]["nseed"])),
            int(float(config["ABC"]["locdim"])),
        ),
    script:
        "../scripts/01.create_abc_simulation_randints.py"


rule create_abc_pod_randints:
    output:
        npy=protected("resources/rand.int.pods.npy"),
    log:
        log1="logs/module01/create_abc_pod_randints.log",
    conda:
        "config/env.yaml"
    params:
        seed=int(float(config["ABC"]["performance"]["pods_seed"])),
        nseed=(
            # number of simulations
            int(float(config["ABC"]["performance"]["nsim"])),
            int(float(config["ABC"]["locdim"])),  # number of loci
            # max number of different pod sets
            int(float(config["ABC"]["performance"]["poddim"]))
        )
    script:
        "../scripts/01.create_abc_pod_randints.py"


rule simulate_treeseqs:
    output:
        tsl=temp("results/abc/simulations/trees/sim_{simid}.locus_{locid}.tsl"),
        params="results/abc/simulations/params/sim_{simid}.locus_{locid}.npy",
    input:
        npy=rules.create_abc_simulation_randints.output.npy,
    log:
        log1="logs/module01/simulate_treeseqs/sim_{simid}.locus_{locid}.log",
    conda:
        "config/env.yaml"
    group:
        "abc_simulation"
    threads: 1
    # resources:
    params:
        mutrate=float(config["ABC"]["simulations"]["mutrate"]),
        recrate=float(config["ABC"]["simulations"]["recrate"]),
        loclen=int(float(config["ABC"]["simulations"]["locuslen"])),
        nsam=int(float(config["ABC"]["simulations"]["nsam"])),
        dtwf=int(float(config["ABC"]["simulations"]["dtwf"])),
        model=config["ABC"]["model"],
    script:
        "../scripts/01.simulate_treeseqs.py"


rule confirm_params_for_loci_being_same_and_save_params:
    output:
        params=protected("results/abc/simulations/params/sim_{simid}.params.npy"),
    input:
        params=expand(
            rules.simulate_treeseqs.output.params,
            locid=wildcards_locid(config),
            allow_missing=True,
        ),
    log:
        log1="logs/module01/confirm_params_for_loci_being_same_and_save_params/sim_{simid}.checked.log",
    conda:
        "config/env.yaml"
    threads: 1
    script:
        "../scripts/01.confirm_params_for_loci_being_same_and_save_params.py"


rule simulate_treeseq_pods:
    output:
        tsl="results/abc/pods/trees/pod_{podid}.tsl",
    input:
        npy=rules.create_abc_pod_randints.output.npy,
    log:
        log1="logs/module01/simulate_treeseq_pods/pod_{podid}.log",
    conda:
        "config/env.yaml"
    threads: 1
    # resources:
    params:
        nsim=int(float(config["ABC"]["performance"]["nsim"])),
        mutrate=float(config["ABC"]["simulations"]["mutrate"]),
        recrate=float(config["ABC"]["simulations"]["recrate"]),
        nloci=int(float(config["ABC"]["simulations"]["nloci"])),
        loclen=int(float(config["ABC"]["simulations"]["locuslen"])),
        nsam=int(float(config["ABC"]["simulations"]["nsam"])),
        dtwf=int(float(config["ABC"]["simulations"]["dtwf"])),
        params=(
            int(float(config["ABC"]["performance"]["pods"][0][0])),
            float(config["ABC"]["performance"]["pods"][1][0]),
            float(config["ABC"]["performance"]["pods"][2][0]),
            int(float(config["ABC"]["performance"]["pods"][3][0])),
        )
    script:
        "../scripts/01.simulate_treeseq_pods.py"


rule generate_discretizing_breakpoints_for_sumstats:
    output:
        npytxt=protected(
            expand(
                "resources/discretization/sumstat_{sumstat}.npytxt",
                sumstat=set(
                    itertools.chain.from_iterable(
                        [
                            sumstat_set.split("/")
                            for sumstat_set in config["ABC"]["sumstats"]
                        ]
                    )
                ),
            )
        ),
    input:
        tsl=expand(rules.simulate_treeseq_pods.output.tsl,
            podid=0)
    log:
        log1="logs/module01/generate_discretizing_breakpoints_for_sumstats.log",
    conda:
        "config/env.yaml"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"]))
    script:
        "../scripts/01.generate_discretizing_breakpoints_for_sumstats.py"


rule calculate_sumstats:
    output:
        npy="results/abc/simulations/sumstats/sim_{simid}.locus_{locid}.npy",
        sumstat_count="results/abc/simulations/sumstats/sim_{simid}.locus_{locid}.names.npy"
    input:
        tsl=rules.simulate_treeseqs.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
    log:
        log1="logs/module01/calculate_sumstats/sim_{simid}.locus_{locid}.log",
    conda:
        "config/env.yaml"
    group:
        "abc_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"]))
    script:
        "../scripts/01.calculate_sumstats.py"


rule calculate_masked_sumstats:
    output:
        npy="results/abc/simulations/sumstats_masked/sim_{simid}.locus_{locid}.npy",
        sumstat_count="results/abc/simulations/sumstats_masked/sim_{simid}.locus_{locid}.names.npy"
    input:
        tsl=rules.simulate_treeseqs.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        mask="resources/mask/locus_{locid}.txt",  # mask file
    log:
        log1="logs/module01/calculate_masked_sumstats/sim_{simid}.locus_{locid}.log",
    conda:
        "config/env.yaml"
    group:
        "abc_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"]))
    script:
        "../scripts/01.calculate_masked_sumstats.py"


rule aggregate_sumstats:
    output:
        sumstats="results/abc/simulations/sumstats.feather",
    input:
        npys=expand(
            rules.calculate_sumstats.output.npy,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        sumstat_names=expand(
            rules.calculate_sumstats.output.sumstat_count,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        params=expand(
            "results/abc/simulations/params/sim_{simid}.params.npy",
            simid=wildcards_simid(config),
        )
    log:
        log1="logs/module01/aggregate_sumstats.log",
    conda:
        "config/env.yaml"
    threads: 1
    # resources:
    params:
        simid_wc = wildcards_simid(config),
        locid_wc = wildcards_locid(config)
    script:
        "../scripts/01.aggregate_sumstats.py"


rule aggregate_sumstats_masked:
    output:
        masked="results/abc/simulations/sumstats.masked.feather"
    input:
        npys_masked=expand(
            "results/abc/simulations/sumstats_masked/sim_{simid}.locus_{locid}.npy",
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        sumstat_names=expand(
            rules.calculate_masked_sumstats.output.sumstat_count,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
            ),
        params=expand(
            "results/abc/simulations/params/sim_{simid}.params.npy",
            simid=wildcards_simid(config),
        ),
    log:
        log1="logs/module01/aggregate_sumstats_masked.log",
    threads: 1
    # resources:
    # params:
    run:
        sys.exit(
            "#" * 600
        + "aggregate_sumstats_masked\n"
            + "here, aggregate parameters and average the sumstates of the loci"
        )


rule alternative_model_create_simulation_randints:
    output:
        npy=protected("resources/rand.alternative_model.int.npy"),
    log:
        log1="logs/module01/alternative_model_create_simulation_randints.log",
    conda:
        "config/env.yaml"
    params:
        seed=int(float(config["ABC"]["alternative_model"]["seed"])),
        nseed=(
            int(float(config["ABC"]["nseed"])),
            int(float(config["ABC"]["locdim"])),
        ),
    script:
        "../scripts/01.create_abc_simulation_randints.py"


rule alternative_model_simulate_treeseqs:
    output:
        tsl=temp("results/abc/simulations/alternative_trees/sim_{simid}.locus_{locid}.tsl"),
        params="results/abc/simulations/alternative_params/sim_{simid}.locus_{locid}.npy",
    input:
        npy=rules.alternative_model_create_simulation_randints.output.npy,
    log:
        log1="logs/module01/alternative_model_simulate_treeseqs/sim_{simid}.locus_{locid}.log",
    conda:
        "config/env.yaml"
    group:
        "alternative_simulation"
    threads: 1
    # resources:
    params:
        mutrate=float(config["ABC"]["simulations"]["mutrate"]),
        recrate=float(config["ABC"]["simulations"]["recrate"]),
        loclen=int(float(config["ABC"]["simulations"]["locuslen"])),
        nsam=int(float(config["ABC"]["simulations"]["nsam"])),
        dtwf=int(float(config["ABC"]["simulations"]["dtwf"])),
        model=config["ABC"]["alternative_model"],
    script:
        "../scripts/01.alternative_model_simulate_treeseqs.py"


rule alternative_model_calculate_sumstats:
    output:
        npy="results/abc/simulations/alternative_sumstats/sim_{simid}.locus_{locid}.npy",
        sumstat_count="results/abc/simulations/alternative_sumstats/sim_{simid}.locus_{locid}.names.npy"
    input:
        tsl=rules.alternative_model_simulate_treeseqs.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
    log:
        log1="logs/module01/alternative_model_calculate_sumstats/sim_{simid}.locus_{locid}.log",
    conda:
        "config/env.yaml"
    group:
        "alternative_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"]))
    script:
        "../scripts/01.calculate_sumstats.py"


rule alternative_model_calculate_masked_sumstats:
    output:
        npy="results/abc/simulations/alternative_sumstats_masked/sim_{simid}.locus_{locid}.npy",
        sumstat_count="results/abc/simulations/alternative_sumstats_masked/sim_{simid}.locus_{locid}.names.npy"
    input:
        tsl=rules.alternative_model_simulate_treeseqs.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        mask="resources/mask/locus_{locid}.txt",  # mask file
    log:
        log1="logs/module01/alternative_model_calculate_masked_sumstats/sim_{simid}.locus_{locid}.log",
    conda:
        "config/env.yaml"
    group:
        "alternative_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"]))
    script:
        "../scripts/01.calculate_masked_sumstats.py"


rule confirm_alternative_params_for_loci_being_same_and_save_params:
    output:
        params=protected("results/abc/simulations/alternative_params/sim_{simid}.params.npy"),
    input:
        params=expand(
            rules.alternative_model_simulate_treeseqs.output.params,
            locid=wildcards_locid(config),
            allow_missing=True,
        ),
    log:
        log1="logs/module01/confirm_alternative_params_for_loci_being_same_and_save_params/sim_{simid}.checked.log",
    conda:
        "config/env.yaml"
    threads: 1
    script:
        "../scripts/01.confirm_params_for_loci_being_same_and_save_params.py"


rule alternative_model_aggregate_sumstats:
    output:
        sumstats="results/abc/simulations/alternative.sumstats.feather",
    input:
        npys=expand(
            rules.alternative_model_calculate_sumstats.output.npy,
            simid=wildcards_simid(config, alternative_model=True),
            locid=wildcards_locid(config),
        ),
        sumstat_names=expand(
            rules.alternative_model_calculate_sumstats.output.sumstat_count,
            simid=wildcards_simid(config, alternative_model=True),
            locid=wildcards_locid(config),
        ),
        params=expand(
            "results/abc/simulations/alternative_params/sim_{simid}.params.npy",
            simid=wildcards_simid(config, alternative_model=True),
        )
    log:
        log1="logs/module01/alternative_model_aggregate_sumstats.log",
    conda:
        "config/env.yaml"
    threads: 1
    # resources:
    params:
        simid_wc = wildcards_simid(config, alternative_model=True),
        locid_wc = wildcards_locid(config)
    script:
        "../scripts/01.aggregate_sumstats.py"


rule alternative_model_aggregate_sumstats_masked:
    output:
        masked="results/abc/simulations/alternative.sumstats.masked.feather"
    input:
        npys_masked=expand(
            rules.alternative_model_calculate_masked_sumstats.output.npy,
            simid=wildcards_simid(config, alternative_model=True),
            locid=wildcards_locid(config),
        ),
        sumstat_names=expand(
            rules.alternative_model_calculate_masked_sumstats.output.sumstat_count,
            simid=wildcards_simid(config, alternative_model=True),
            locid=wildcards_locid(config),
            ),
        params=expand(
            "results/abc/simulations/alternative_params/sim_{simid}.params.npy",
            simid=wildcards_simid(config, alternative_model=True),
        ),
    log:
        log1="logs/module01/alternative_model_aggregate_sumstats_masked.log",
    threads: 1
    # resources:
    # params:
    run:
        sys.exit(
            "#" * 600
        + "alternative_model_aggregate_sumstats_masked\n"
            + "here, aggregate parameters and average the sumstates of the loci"
        )


rule calculate_podstats:
    output:
        npy="results/abc/pods/podstats/podid_{podid}.npy",
        sumstat_count="results/abc/pods/podstats/podid_{podid}.names.npy"
    input:
        tsl=rules.simulate_treeseq_pods.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
    log:
        log1="logs/module01/calculate_podstats/podid_{podid}.log",
    conda:
        "config/env.yaml"
    group:
        "abc_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"]))
    script:
        "../scripts/01.calculate_podstats.py"


rule aggregate_podstats:
    output:
        sumstats="results/abc/pods/podstats.feather",
    input:
        npys=expand(rules.calculate_podstats.output.npy,
            podid=wildcards_podid(config)),
        sumstat_names=expand(
            rules.calculate_podstats.output.sumstat_count,
            podid=wildcards_podid(config)),
    log:
        log1="logs/module01/aggregate_podstats.log",
    conda:
        "config/env.yaml"
    threads: 1
    # resources:
    params:
        podid_wc=wildcards_podid(config),
    script:
        "../scripts/01.aggregate_podstats.py"


rule calculate_podstats_masked:
    output:
        npy="results/abc/pods/podstats_masked/podid_{podid}.npy",
        sumstat_count="results/abc/pods/podstats_masked/podid_{podid}.names.npy"
    input:
        tsl=rules.simulate_treeseq_pods.output.tsl,


rule aggregate_podstats_masked:
    output:
        sumstats="results/abc/pods/podstats.masked.feather",
    input:
        sumstats=expand(rules.calculate_podstats_masked.output.npy,
            podid=wildcards_podid(config)),
        sumstat_names=expand(
            rules.calculate_podstats_masked.output.sumstat_count,
            podid=wildcards_podid(config)),
    log:
        log1="logs/module01/aggregate_podstats_masked.log",
    #conda:
    #    "config/env.yaml"
    threads: 1
    # resources:
    params:
        podid_wc=wildcards_podid(config),
    run:
        sys.exit(
            "#" * 600
        + "aggregate_podstats\n"
            + "here, aggregate parameters and average the sumstates of the loci"
        )

