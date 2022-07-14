"""
The helper module of the workflow

Here, we create some rules, that are not needed for the main workflow, but
provided useful and needed information to start proper and meaningful running
of the workflow
"""


rule calculate_theta_watterson:
    output:
        txt="theta_watterson_{mode}",
    input:
        tsl_pod=expand("results/abc/pods/trees/pod_{podid}.tsl", podid=0),
        treeseq_athal=config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["path"],
        samples=config["ABC"]["athaliana"]["observations"]["treeseq_1001"]["popid"][
            config["ABC"]["sumstats_specs"]["TM_WIN"]["which_pop_if_athal"]
        ]["path"],
    log:
        log1="logs/module00/calculate_two_N_zero.{mode}.log",
    conda:
        "../../config/env.yaml"
    params:
        seed=int(
            float(
                config["ABC"]["sumstats_specs"]["TM_WIN"][
                    "two_N_zero_if_expected_seed"
                ]
            )
        ),
        chrom_multiplier=float(
            config["ABC"]["athaliana"]["observations"]["treeseq_1001"][
                "chrom_multiplier"
            ]
        ),
    script:
        "../scripts/00.calculate_theta_watterson.py"
