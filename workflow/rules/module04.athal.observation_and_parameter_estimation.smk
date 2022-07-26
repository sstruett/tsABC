"""
Here, we provide a pipeline which runs the whole abc on Arabidopsis thaliana
data. We observe the statistics from Arabidopsis thaliana same as in the abc
module. It makes sense to provide the first POD in the first model to be
somewhat inspired to the values we hopefulle know already for thaliana. In a
second part, we run simulations under a six-parameter-model and perform the
parameter estimation on it.
"""


localrules:
    confirm_athal_params_for_loci_being_same_and_save_params,


rule create_athal_simulation_randints:
    output:
        npy=protected("resources/rand.int.athal.npy"),
    log:
        log1="logs/module04/create_athal_simulation_randints.log",
    conda:
        "../config/env.yaml"
    params:
        seed=int(float(config["ABC"]["athaliana"]["seed"])),
        nseed=(
            int(float(config["ABC"]["athaliana"]["nseed"])),
            int(float(config["ABC"]["athaliana"]["locdim"])),
        ),
    script:
        "../scripts/01.create_abc_simulation_randints.py"


rule simulate_athal_treeseqs:
    output:
        tsl=temp("results/athal/simulations/trees/sim_{simid}.locus_{locid}.tsl"),
        params="results/athal/simulations/params/sim_{simid}.locus_{locid}.npy",
    input:
        npy=rules.create_athal_simulation_randints.output.npy,
    log:
        log1="logs/module04/simulate_athal_treeseqs/sim_{simid}.locus_{locid}.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_simulation"
    threads: 1
    # resources:
    params:
        mutrate=float(config["ABC"]["simulations"]["mutrate"]),
        recrate=float(config["ABC"]["simulations"]["recrate"]),
        loclen=int(float(config["ABC"]["simulations"]["locuslen"])),
        nsam=int(float(config["ABC"]["simulations"]["nsam"])),
        dtwf=int(float(config["ABC"]["simulations"]["dtwf"])),
        model=config["ABC"]["athaliana"],
    script:
        "../scripts/04.simulate_sixparmodel_treeseqs.py"


rule confirm_athal_params_for_loci_being_same_and_save_params:
    output:
        params=protected("results/athal/simulations/params/sim_{simid}.params.npy"),
    input:
        params=expand(
            rules.simulate_athal_treeseqs.output.params,
            locid=wildcards_locid(config),
            allow_missing=True,
        ),
    log:
        log1="logs/module04/confirm_athal_params_for_loci_being_same_and_save_params/sim_{simid}.checked.log",
    conda:
        "../../config/env.yaml"
    threads: 1
    script:
        "../scripts/01.confirm_params_for_loci_being_same_and_save_params.py"


rule calculate_athal_sumstats:
    output:
        npy="results/athal/simulations/sumstats/sim_{simid}.locus_{locid}.npy",
        sumstat_count="results/athal/simulations/sumstats/sim_{simid}.locus_{locid}.names.npy",
    input:
        tsl=rules.simulate_athal_treeseqs.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
    log:
        log1="logs/module04/calculate_athal_sumstats/sim_{simid}.locus_{locid}.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed_athal"])),
    script:
        "../scripts/01.calculate_sumstats.py"


rule calculate_athal_masked_sumstats:
    output:
        npy="results/athal/simulations/sumstats_masked/sim_{simid}.locus_{locid}.npy",
        sumstat_count="results/athal/simulations/sumstats_masked/sim_{simid}.locus_{locid}.names.npy",
    input:
        tsl=rules.simulate_athal_treeseqs.output.tsl,
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        mask="resources/mask/locus_{locid}.txt",  # mask file
    log:
        log1="logs/module04/calculate_athal_masked_sumstats/sim_{simid}.locus_{locid}.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_simulation"
    threads: 1
    # resources:
    params:
        sumstats=config["ABC"]["sumstats"],
        seed=int(float(config["ABC"]["sumstats_specs"]["seed"])),
    script:
        "../scripts/01.calculate_masked_sumstats.py"


rule aggregate_athal_sumstats:
    output:
        sumstats="results/athal/simulations/sumstats.feather",
    input:
        npys=expand(
            rules.calculate_athal_sumstats.output.npy,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        sumstat_names=expand(
            rules.calculate_athal_sumstats.output.sumstat_count,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        params=expand(
            "results/athal/simulations/params/sim_{simid}.params.npy",
            simid=wildcards_simid(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_sumstats.log",
    conda:
        "../../config/env.yaml"
    threads: 1
    # resources:
    params:
        simid_wc=wildcards_simid(config),
        locid_wc=wildcards_locid(config),
    script:
        "../scripts/01.aggregate_sumstats.py"


rule aggregate_athal_sumstats_masked:
    output:
        sumstats="results/athal/simulations/sumstats.masked.feather",
    input:
        npys=expand(
            rules.calculate_athal_masked_sumstats.output.npy,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        sumstat_names=expand(
            rules.calculate_masked_sumstats.output.sumstat_count,
            simid=wildcards_simid(config),
            locid=wildcards_locid(config),
        ),
        params=expand(
            rules.confirm_athal_params_for_loci_being_same_and_save_params.output.params,
            simid=wildcards_simid(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_sumstats_masked.log",
    conda:
        "../../config/env.yaml"
    threads: 1
    # resources:
    params:
        simid_wc = wildcards_simid(config),
        locid_wc = wildcards_locid(config)
    script:
        "../scripts/01.aggregate_sumstats.py"


rule calculate_athal_observations:
    output:
        sumstats="results/athal/observed/population_{popid}.unmasked.feather",
    input:
        athal_treeseq=config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["path"],
        sample_list=lambda wildcards: config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["popid"][wildcards.popid]["path"],
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        maskfiles=expand("resources/mask/locus_{chromid}.txt",
            chromid=range(5),  # the chromosome ids are 0 based
            )
    log:
        log1="logs/module04/calculate_athal_observations/population_{popid}.log",
    conda:
        "../../config/env.yaml"
    threads: 1
    params:
        masked=False,
        chrom_multiplier=float(config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["chrom_multiplier"]),
        seed=int(float(config["ABC"]["athaliana"]["observations"]["seed"])),
        sumstats=config["ABC"]["sumstats"],
    script:
        "../scripts/04.calculate_athal_observations.py"


rule calculate_athal_observations_masked:
    output:
        sumstats="results/athal/observed/population_{popid}.masked.feather"
    input:
        athal_treeseq=config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["path"],
        sample_list=lambda wildcards: config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["popid"][wildcards.popid]["path"],
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        maskfiles=expand("resources/mask/locus_{chromid}.txt",
            chromid=range(5),  # the chromosome ids are 0 based
            )
    log:
        log1="logs/module04/calculate_athal_observations_masked/population_{popid}.log"
    conda:
        "../../config/env.yaml"
    threads: 1
    params:
        masked=True,
        chrom_multiplier=float(config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["chrom_multiplier"]),
        seed=int(float(config["ABC"]["athaliana"]["observations"]["seed"])),
        sumstats=config["ABC"]["sumstats"],
    script:
        "../scripts/04.calculate_athal_observations.py"


rule calculate_athal_all_region_observations:
    output:
        sumstats="results/athal/observed/population_{popid}.all_region.feather"
    input:
        athal_treeseq=config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["path"],
        sample_list=lambda wildcards: config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["popid"][wildcards.popid]["path"],
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        # maskfiles must be provided in chromosome sequence
        maskfiles=expand("resources/mask/locus_{chromid}.txt",
            chromid=range(5),  # the chromosome ids are 0 based
            )
    log:
        log1="logs/module04/calculate_athal_all_region_observations/population_{popid}.log"
    conda:
        "../../config/env.yaml"
    threads: 1
    params:
        masked=False,
        chrom_multiplier=float(config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["chrom_multiplier"]),
        seed=int(float(config["ABC"]["athaliana"]["observations"]["seed"])),
        sumstats=config["ABC"]["sumstats"],
    script:
        "../scripts/04.calculate_athal_all_region_observations.py"


rule calculate_athal_all_region_observations_masked:
    output:
        sumstats="results/athal/observed/population_{popid}.all_region.masked.feather"
    input:
        athal_treeseq=config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["path"],
        sample_list=lambda wildcards: config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["popid"][wildcards.popid]["path"],
        breakpoints=rules.generate_discretizing_breakpoints_for_sumstats.output.npytxt,
        # maskfiles must be provided in chromosome sequence
        maskfiles=expand("resources/mask/locus_{chromid}.txt",
            chromid=range(5),  # the chromosome ids are 0 based
            )
    log:
        log1="logs/module04/calculate_athal_all_region_observations_masked/population_{popid}.log"
    conda:
        "../../config/env.yaml"
    threads: 1
    params:
        masked=True,
        chrom_multiplier=float(config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["chrom_multiplier"]),
        seed=int(float(config["ABC"]["athaliana"]["observations"]["seed"])),
        sumstats=config["ABC"]["sumstats"],
    script:
        "../scripts/04.calculate_athal_all_region_observations.py"


rule aggregate_athal_observations:
    output:
        sumstats="results/athal/observed/sumstats.feather",
        # maps the sumstats (row) to all the specifications of the observation
        identifier="results/athal/observed/dataset.map.feather",
    input:
        observations=expand(
            "results/athal/observed/population_{popid}.unmasked.feather",
            popid=wildcards_popid(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_observations.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_observation"
    threads: 1
    script:
        "../scripts/04.aggregate_athal_sumstats.py"


rule aggregate_athal_observations_masked:
    output:
        sumstats="results/athal/observed/sumstats.masked.feather",
        # maps the sumstats (row) to all the specifications of the observation
        identifier="results/athal/observed/dataset.masked.map.feather",
    input:
        observations=expand(
            "results/athal/observed/population_{popid}.masked.feather",
            popid=wildcards_popid(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_observations_masked.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_observation"
    threads: 1
    script:
        "../scripts/04.aggregate_athal_sumstats.py"


rule aggregate_athal_all_regions_observations:
    output:
        sumstats="results/athal/observed/sumstats.all_regions.feather",
        # maps the sumstats (row) to all the specifications of the observation
        identifier="results/athal/observed/dataset.all_regions.map.feather",
    input:
        observations=expand(
            "results/athal/observed/population_{popid}.all_region.feather",
            popid=wildcards_popid(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_all_regions_observations.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_observation"
    threads: 1
    script:
        "../scripts/04.aggregate_athal_sumstats.py"


rule aggregate_athal_all_regions_observations_masked:
    output:
        sumstats="results/athal/observed/sumstats.all_regions.masked.feather",
        # maps the sumstats (row) to all the specifications of the observation
        identifier="results/athal/observed/dataset.all_regions.masked.map.feather",
    input:
        observations=expand(
            "results/athal/observed/population_{popid}.all_region.masked.feather",
            popid=wildcards_popid(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_all_regions_observations_masked.log",
    conda:
        "../../config/env.yaml"
    group:
        "athal_observation"
    threads: 1
    script:
        "../scripts/04.aggregate_athal_sumstats.py"


rule athal_parameter_estimation:
    output:
        estims="results/athal/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}..RDS",
        postplots=directory(
            "results/athal/parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}/"
        ),
    input:
        sumstats="results/athal/transformation/statcomp_{statcomposition}..simulations_sumstats.txt",
        observed="results/athal/transformation/statcomp_{statcomposition}..observed_sumstats.txt",
        identifier=rules.aggregate_athal_observations.output.identifier,
    log:
        log1="logs/module04/athal_parameter_estimation/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "../../config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
        regression="loclinear",  # r-abc package
    script:
        "../scripts/04.parameter_estimation.R"


rule athal_parameter_estimation_masked:
    output:
        estims="results/athal/parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}...RDS",
        postplots=directory(
            "results/athal/parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}/"
        ),
    input:
        sumstats="results/athal/transformation/statcomp_{statcomposition}..simulations_sumstats.masked.txt",
        observed="results/athal/transformation/statcomp_{statcomposition}..observed_sumstats.masked.txt",
        identifier=rules.aggregate_athal_observations_masked.output.identifier,
    log:
        log1="logs/module04/athal_parameter_estimation_masked/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "../../config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
        regression="loclinear",  # r-abc package
    script:
        "../scripts/04.parameter_estimation.R"


rule athal_parameter_estimation_all_regions:
    output:
        estims="results/athal/parameter_estimation_allregions/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}..RDS",
        postplots=directory(
            "results/athal/parameter_estimation_allregions/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}/"
        ),
    input:
        sumstats="results/athal/transformation/statcomp_{statcomposition}..simulations_sumstats.txt",
        observed="results/athal/transformation/statcomp_{statcomposition}..observed_sumstats.all_regions.txt",
        identifier=rules.aggregate_athal_all_regions_observations.output.identifier,
    log:
        log1="logs/module04/athal_parameter_estimation_all_regions/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "../../config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
        regression="loclinear",  # r-abc package
    script:
        "../scripts/04.parameter_estimation.R"


rule athal_parameter_estimation_masked_all_regions:
    output:
        estims="results/athal/parameter_estimation_masked_all_regions/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}...RDS",
        postplots=directory(
            "results/athal/parameter_estimation_masked_all_regions/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}/"
        ),
    input:
        sumstats="results/athal/transformation/statcomp_{statcomposition}..simulations_sumstats.masked.txt",
        observed="results/athal/transformation/statcomp_{statcomposition}..observed_sumstats.all_regions.masked.txt",
        identifier=rules.aggregate_athal_all_regions_observations_masked.output.identifier,
    log:
        log1="logs/module04/athal_parameter_estimation_masked_all_regions/statcomp_{statcomposition}.pls_{plsid}.tolid_{tolid}.transformed.log",
    conda:
        "../../config/env.yaml"
    params:
        ntol=lambda wildcards: int(float(wildcards.tolid)),
        pls=lambda wildcards: int(float(wildcards.plsid)),
        regression="loclinear",  # r-abc package
    script:
        "../scripts/04.parameter_estimation.R"


rule aggregate_athal_parameter_estimation:
    output:
        estims="results/athal/parameter_estimation.RDS",
    input:
        estims=expand(
            rules.athal_parameter_estimation.output.estims,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ) + expand(
            rules.athal_parameter_estimation_all_regions.output.estims,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_parameter_estimation.log",
    # conda:
    #    "../../config/env.yaml"
    # params:
    script:
        "../scripts/04.aggregate_athal_parameter_estimation.R"


rule aggregate_athal_parameter_estimation_masked:
    output:
        estims="results/athal/parameter_estimation_masked.RDS",
    input:
        estims=expand(
            rules.athal_parameter_estimation_masked.output.estims,
            plsid=wildcards_plsid(config),
            tolid=wildcards_tolid(config),
            statcomposition=wildcards_statcomposition(config),
        ),
    log:
        log1="logs/module04/aggregate_athal_parameter_estimation_masked.log",
    # conda:
    #    "../../config/env.yaml"
    # params:
    script:
        "../scripts/04.aggregate_athal_parameter_estimation.R"


rule summarize_athal_param_estims:
    output:
        table="results/athal/summarized_athal_parameter_estimations.csv",
        plot="results/athal/summarized_athal_parameter_estimations.pdf",
    input:
        estims=rules.aggregate_athal_parameter_estimation.output.estims,
    log:
        log1="logs/module04/summarize_athal_param_estims.log",
    # conda:
    #    "../../config/env.yaml"
    # params:
    run:
        sys.exit(
            "#" * 600 + "inside aggregate_athal_parameter_estimation_masked\n" + ""
        )


rule summarize_athal_param_estims_masked:
    output:
        table="results/athal/summarized_athal_parameter_estimations_masked.csv",
        plot="results/athal/summarized_athal_parameter_estimations_masked.pdf",
    input:
        estims=rules.aggregate_athal_parameter_estimation_masked.output.estims,
    log:
        log1="logs/module04/summarize_athal_param_estims_masked.log",
    # conda:
    #    "../../config/env.yaml"
    # params:
    run:
        sys.exit(
            "#" * 600 + "inside aggregate_athal_parameter_estimation_masked\n" + ""
        )


rule transform_athal_sumstats:
    output:
        sumstats="results/athal/transformation/statcomp_{statcomposition}..{simorpod}_{dataset}.txt",
    input:
        sumstats="results/athal/{simorpod}/{dataset}.feather",
        transformer=which_transformer_athal,
        script="workflow/scripts/transformer",
        subsetter="workflow/scripts/subset_table_and_separate_params.py",
    log:
        log1="logs/module04/transform_athal_sumstats/transformation_{statcomposition}.{simorpod}.{dataset}.log",
    conda:
        "../../config/env.yaml"
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


rule find_athal_pls:
    output:
        transformer="results/athal/transformation/transformer_{statcomposition}/Routput_{statcomposition}",
        rmse_plot="results/athal/transformation/transformer_{statcomposition}/RMSE_{statcomposition}.pdf",
    input:
        sumstats=rules.aggregate_athal_sumstats.output.sumstats,
        script="workflow/scripts/find_pls_by_wanted_sumstats.R",
    log:
        log1="logs/module04/find_athal_pls/transformed_{statcomposition}.log",
    conda:
        "../../config/env.yaml"
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


rule find_athal_pls_masked:
    output:
        transformer="results/athal/transformation/masked_transformer_{statcomposition}/Routput_{statcomposition}",
        rmse_plot="results/athal/transformation/masked_transformer_{statcomposition}/RMSE_{statcomposition}.pdf",
    input:
        sumstats=rules.aggregate_athal_sumstats_masked.output.sumstats,
        script="workflow/scripts/find_pls_by_wanted_sumstats.R",
    log:
        log1="logs/module04/find_athal_pls_masked/transformed_{statcomposition}.log",
    conda:
        "../../config/env.yaml"
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
