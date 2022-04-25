"""
The helper module of the workflow

Here, we create the target file list and print them to the standard error
stream. So, we still maintain creating the dag/rulegraph by piping to dot
correct.

Additionally, we intend this module to contain all helper functions which are
needed for the workflow in order to maintain readability in the parent
Snakefile.
"""

import sys
import warnings
import math


def target_files(wildcards, verbose=True):
    """
    function to provide the target file as a list dependent on the wildcards
    """
    target_file_list = []

    # aggregated abc simulations, sumstats and masked sumstats
    if True:
        target_file_list.append("results/abc/simulations/sumstats.feather")

        if config["ABC"]["sumstats_specs"]["run_masked"]:
            target_file_list.append("results/abc/simulations/sumstats.masked.feather")

    # only simulate alternative model if it asked for
    if config["ABC"]["alternative_model"]["do_model_choice"]:
        # simulations
        target_file_list.append("results/abc/simulations/alternative.sumstats.feather")

        # model choice results
        target_file_list.append("results/abc/model_choice.RDS")

        if config["ABC"]["sumstats_specs"]["run_masked"]:
            target_file_list.append(
                "results/abc/simulations/alternative.sumstats.masked.feather"
            )


    # parameter estimations
    if True:
        target_file_list.append("results/abc/parameter_estimation.RDS")

        # if there is the masking asked for
        if config["ABC"]["sumstats_specs"]["run_masked"]:
            target_file_list.append("results/abc/parameter_estimation_masked.RDS")


    # print the requested files to the standard error stream
    if verbose:
        print("\n" + "_" * 80, file=sys.stderr)
        for file_index, target_file in enumerate(target_file_list, start=1):
            print(f"  {file_index}.) {target_file}", file=sys.stderr)
        print("\n" + "=" * 80 + "\n" * 2, file=sys.stderr)

    return target_file_list


def check_configuration_file(config):
    """
    check a few necessities on some values provided in the config.yaml
    """
    assert float(config["ABC"]["nseed"]) >= float(
        config["ABC"]["simulations"]["nsim"]
    ), "you did not create enough seeds, error 1"
    assert float(config["ABC"]["locdim"]) >= float(
        config["ABC"]["simulations"]["nloci"]
    ), "you did not create enough seeds, error 2"

    # check prior definitions

    # loop through priors of model
    for prior in config["ABC"]["model"]["priors"]:
        if prior[3] == "loguniform":
            assert float(prior[0]) > 0, "".join(
                [
                    "technically possible, but a prior range low of 0 or smaller",
                    " provides some problems that have to be checked and are ",
                    "currently not implemented",
                ]
            )
        assert float(prior[0]) <= float(prior[1]), "maldefined prior range"
        assert prior[2] in ["int", "float",], "".join(
            [
                "maldefined prior numeric type, only 'int'",
                " or 'float' are supported so far",
            ]
        )
        assert prior[3] in ["uniform", "loguniform",], "".join(
            [
                "maldefined prior distribution type; ",
                "so far supported: 'uniform' or 'loguniform'",
            ]
        )

    # loop through priors of alternative model
    for prior in config["ABC"]["alternative_model"]["priors"]:
        if prior[3] == "loguniform":
            assert float(prior[0]) > 0, "".join(
                [
                    "technically possible, but a prior range low of 0 or smaller",
                    " provides some problems that have to be checked and are ",
                    "currently not implemented",
                ]
            )
        assert float(prior[0]) <= float(prior[1]), "maldefined prior range"
        assert prior[2] in ["int", "float",], "".join(
            [
                "maldefined prior numeric type, only 'int'",
                " or 'float' are supported so far",
            ]
        )
        assert prior[3] in ["uniform", "loguniform",], "".join(
            [
                "maldefined prior distribution type; ",
                "so far supported: 'uniform' or 'loguniform'",
            ]
        )

    # check that all pod definitions are complete (whether pods params are of
    # same length)
    pod_param_list_lenghts = [
        len(pod_param_list) for pod_param_list in config["ABC"]["performance"]["pods"]
    ]
    assert (
        len(list(set(pod_param_list_lenghts))) == 1
    ), "pod parameters do not have the same length"

    # sumstats and sumstat sets
    listed_sumstats = set()
    for sumstat_set in config["ABC"]["sumstats"]:
        listed_sumstats.update(sumstat_set.split("/"))
    for singel_sumstat in listed_sumstats:
        assert singel_sumstat in {"SFS", "TM_WIN", "LD",}, "".join(
            [
                "the calculation of the summarizing",
                f" statistic '{singel_sumstat}' has not yet been implemented",
            ]
        )

    # check that run_masked in config file is bool
    assert (
        type(config["ABC"]["sumstats_specs"]["run_masked"]) == bool
    ), 'wrong type in config file config["ABC"]["sumstats_specs"]["run_masked"] must be boolean'

    # check that do_model_choice in config file is bool
    assert (
        type(config["ABC"]["alternative_model"]["do_model_choice"]) == bool
    ), 'wrong type in config file config["ABC"]["sumstats_specs"]["run_masked"] must be boolean'

    return True


def wildcards_simid(config, alternative_model=False):
    """
    provide a list/generator of all simulation ids. This a helper function for
    the aggregation rules; it solely relies on the configuration yaml. Here the
    sim ids may stand for a set of simulated, in case in the configuration file
    jobcluster_simulations is set

    Args:
        config: dict; configuration dict, usually provided from snakemake which
            reads in as yaml
        alternative_model: bool, whether to get the simids for the alternative
            model simulations or not; default False

    Returns:
        list of simids (int)
    """
    if alternative_model:
        nsim = int(float(config["ABC"]["alternative_model"]["nsim"]))
        low = int(float(config["ABC"]["alternative_model"]["first_simid"]))
        high = low + nsim

        # number of aggregated sims, i. e. the cluster size
        jclus = int(float(config["ABC"]["alternative_model"]["jobluster_simulations"]))

        # set the simids named as the first for each clustered set
        simids = [low + clusid * jclus for clusid in range(math.ceil(nsim / jclus))]
    else:
        # range of the actual simids
        nsim = int(float(config["ABC"]["simulations"]["nsim"]))
        low = int(float(config["ABC"]["first_simid"]))
        high = low + nsim

        # number of aggregated sims, i. e. the cluster size
        jclus = int(float(config["ABC"]["jobluster_simulations"]))

        # set the simids named as the first for each clustered set
        simids = [low + clusid * jclus for clusid in range(math.ceil(nsim / jclus))]

    return simids


def wildcards_locid(config):
    """
    provide a list/generator of all locus ids. This a helper function for
    the aggregation rules; it solely relies on the configuration yaml.
    """
    return range(int(float(config["ABC"]["simulations"]["nloci"])))


def wildcards_podid(config):
    """
    provide a list/generator of all pseudo-observed dataset ids. This a helper
    function for the aggregation rules; it solely relies on the configuration
    yaml.
    """
    return range(len(config["ABC"]["performance"]["pods"]))


def wildcards_plsid(config):
    """
    provide a list/generator of all pls ids. This a helper function for the
    aggregation rules; it solely relies on the configuration yaml.
    """
    return [int(float(plsid)) for plsid in config["ABC"]["performance"]["pls"]]


def wildcards_tolid(config):
    """
    provide a list/generator of all tolerance ids. This a helper function for
    the aggregation rules; it solely relies on the configuration yaml.
    """
    return [
        int(float(tolid))
        for tolid in config["ABC"]["performance"]["tolerance_estimates"]
    ]


def wildcards_statcomposition(config):
    """
    provide a list/generator of all statistics composition ids. This a helper
    function for the aggregation rules; it solely relies on the configuration
    yaml.
    """
    return [
        statcompositionid
        for statcompositionid, _ in enumerate(config["ABC"]["sumstats"])
    ]


def which_transformer(wildcards):
    """
    based on the wildcards of the rule, the correct transformer must be found;
    this basically means if it is the transformer for the masked file or not
    """
    if any(["masked" in wildcard for wildcard in wildcards]):
        transformer = f"results/abc/transformation/masked_transformer_{wildcards['statcomposition']}/Routput_{wildcards['statcomposition']}"
    else:
        transformer = f"results/abc/transformation/transformer_{wildcards['statcomposition']}/Routput_{wildcards['statcomposition']}"

    return transformer