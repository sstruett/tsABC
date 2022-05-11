"""
Functions for the ABC
"""


import sys
import datetime
import numpy as np
import tskit
import msprime
import itertools
import warnings
import math


def simulate_treesequence_under_model(params, rule_parameters, rng, log):
    """Simulates a tree sequence

    The simulation model is defined in this function. The simulation takes place
    using the provided parameters. We only need two phases for the change in
    effective recombination rate. The change of a model is implemented already.

    Args:
        params: A list containing the provided parameters. Here we use a 4
            parameter model
        rule_parameters: Non random parameters of the model
        rng: The provided random number generator
        log: String, filepath to log file, will append some logs to it

    Returns:
        A single tree sequnce, with mutations
    """
    # create demography for phase zero
    demography_phase_zero = msprime.Demography()
    demography_phase_zero.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[1], make_int=True
        ),
        description="Single population",
    )

    # create simulation model for phase zero
    model_phase_zero = [
        msprime.DiscreteTimeWrightFisher(
            duration=min(rule_parameters["dtwf"], params[3])
        ),
        msprime.SmcPrimeApproxCoalescent(),
    ]

    # prepare the paramters for the simulation; this is the model implementation
    param_dict_list = [{}, {}]  # one dict per phase
    # phase 0 - selfing
    param_dict_list[0]["samples"] = rule_parameters["nsam"]
    param_dict_list[0]["demography"] = demography_phase_zero
    param_dict_list[0]["sequence_length"] = int(float(rule_parameters["loclen"]))
    param_dict_list[0][
        "recombination_rate_effective"
    ] = rescale_recombination_rate_by_selfing(
        recombination_rate=rule_parameters["recrate"], selfing_rate=params[1]
    )
    param_dict_list[0]["msprime_model"] = model_phase_zero
    param_dict_list[0]["end_time"] = params[3]

    # create demography for phase one
    demography_phase_one = msprime.Demography()
    demography_phase_one.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[2], make_int=True
        ),
        description="Single population",
    )

    # create simulation model for phase one
    model_phase_one = [
        msprime.DiscreteTimeWrightFisher(
            duration=max(0, rule_parameters["dtwf"] - params[3])
        ),
        msprime.SmcPrimeApproxCoalescent(),
    ]

    # phase 1 - outcrossing
    param_dict_list[1]["demography"] = demography_phase_one
    param_dict_list[1][
        "recombination_rate_effective"
    ] = rescale_recombination_rate_by_selfing(
        recombination_rate=rule_parameters["recrate"], selfing_rate=params[2]
    )
    param_dict_list[1]["msprime_model"] = model_phase_one
    param_dict_list[1]["start_time"] = params[3]

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("created model", file=logfile)

    demography_total = msprime.Demography()
    demography_total.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[1], make_int=True
        ),
        description="Single population",
    )
    demography_total.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[2], make_int=True
        ),
        description="Single population",
    )

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(
            "".join(
                [
                    "\n",
                    "_" * 80,
                    "\nDemography for simulation\n",
                    "This demography is copy pasted from the creation of the",
                    "two phases, which we need to change the recombination rate",
                    "through time.\n\n",
                    str(demography_total),
                    "\n" + "=" * 80,
                ]
            ),
            file=logfile,
        )

    ts = simulate_transition_to_selfing(param_dict_list, rng)

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("simulated ancestry", file=logfile)

    seed = rng.integers(low=1, high=2**32 - 1, size=1)
    ts_mutated = msprime.sim_mutations(
        ts,
        rate=rule_parameters["mutrate"],
        random_seed=seed,
        model=None,
        start_time=None,
        end_time=None,
        discrete_genome=True,
        keep=True,
    )

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("simulated mutations", file=logfile)

    return ts_mutated


def simulate_treesequence_under_alternative_model(params, rule_parameters, rng, log):
    """Simulates a tree sequence

    The simulation model is defined in this function. The simulation takes place
    using the provided parameters. We simulate a single population with constant
    selfing rate undergoing a single stepwise change in population size.

    Args:
        params: A list containing the provided parameters. Here we use a 4
            parameter model
        rule_parameters: Non random parameters of the model
        rng: The provided random number generator
        log: String, filepath to log file, will append some logs to it

    Returns:
        A single tree sequnce, with mutations
    """
    # create demography for phase zero
    demography = msprime.Demography()
    demography.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[1], make_int=True
        ),
        description="Single population",
    )
    demography.add_population_parameters_change(
        params[3],
        initial_size=rescale_population_size_by_selfing(
            population_size=params[2], selfing_rate=params[1], make_int=True
        ),
        growth_rate=None,
        population=None,
    )

    # create simulation model for phase zero
    msprime_model = [
        msprime.DiscreteTimeWrightFisher(
            duration=rule_parameters["dtwf"],
        ),
        msprime.SmcPrimeApproxCoalescent(),
    ]

    # prepare the paramters for the simulation; this is the model implementation
    param_dict = {}
    param_dict["samples"] = rule_parameters["nsam"]
    param_dict["demography"] = demography
    param_dict["sequence_length"] = int(float(rule_parameters["loclen"]))
    param_dict["recombination_rate_effective"] = rescale_recombination_rate_by_selfing(
        recombination_rate=rule_parameters["recrate"], selfing_rate=params[1]
    )
    param_dict["msprime_model"] = msprime_model

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("created model", file=logfile)

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(
            "".join(
                [
                    "\n",
                    "_" * 80,
                    "\nDemography for simulation\n\n",
                    str(demography),
                    "\n" + "=" * 80,
                ]
            ),
            file=logfile,
        )

    ts = simulate_popsize_change_with_constant_selfing(param_dict, rng)

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("simulated ancestry", file=logfile)

    seed = rng.integers(low=1, high=2**32 - 1, size=1)
    ts_mutated = msprime.sim_mutations(
        ts,
        rate=rule_parameters["mutrate"],
        random_seed=seed,
        model=None,
        start_time=None,
        end_time=None,
        discrete_genome=True,
        keep=True,
    )

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("simulated mutations", file=logfile)

    return ts_mutated


def simulate_treesequence_under_six_parameter_model(params, rule_parameters, rng, log):
    """Simulates a tree sequence

    The simulation model is defined in this function. The simulation takes place
    using the provided parameters. We only need two phases for the change in
    effective recombination rate. The change of a model is implemented already.

    Args:
        params: A list containing the provided parameters. Here we use a 6
            parameter model
        rule_parameters: Non random parameters of the model
        rng: The provided random number generator
        log: String, filepath to log file, will append some logs to it

    Returns:
        A single tree sequnce, with mutations
    """
    # create demography for phase zero
    demography_phase_zero = msprime.Demography()
    demography_phase_zero.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[1], make_int=True
        ),
        description="Single population",
    )

    # create simulation model for phase zero
    model_phase_zero = [
        msprime.DiscreteTimeWrightFisher(
            duration=min(rule_parameters["dtwf"], params[3])
        ),
        msprime.SmcPrimeApproxCoalescent(),
    ]

    # prepare the paramters for the simulation; this is the model implementation
    param_dict_list = [{}, {}]  # one dict per phase
    # phase 0 - selfing
    param_dict_list[0]["samples"] = rule_parameters["nsam"]
    param_dict_list[0]["demography"] = demography_phase_zero
    param_dict_list[0]["sequence_length"] = int(float(rule_parameters["loclen"]))
    param_dict_list[0][
        "recombination_rate_effective"
    ] = rescale_recombination_rate_by_selfing(
        recombination_rate=rule_parameters["recrate"], selfing_rate=params[1]
    )
    param_dict_list[0]["msprime_model"] = model_phase_zero
    param_dict_list[0]["end_time"] = params[3]

    # create demography for phase one
    demography_phase_one = msprime.Demography()
    demography_phase_one.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[2], make_int=True
        ),
        description="Single population",
    )

    # create simulation model for phase one
    model_phase_one = [
        msprime.DiscreteTimeWrightFisher(
            duration=max(0, rule_parameters["dtwf"] - params[3])
        ),
        msprime.SmcPrimeApproxCoalescent(),
    ]

    # phase 1 - outcrossing
    param_dict_list[1]["demography"] = demography_phase_one
    param_dict_list[1][
        "recombination_rate_effective"
    ] = rescale_recombination_rate_by_selfing(
        recombination_rate=rule_parameters["recrate"], selfing_rate=params[2]
    )
    param_dict_list[1]["msprime_model"] = model_phase_one
    param_dict_list[1]["start_time"] = params[3]

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("created model", file=logfile)

    demography_total = msprime.Demography()
    demography_total.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[1], make_int=True
        ),
        description="Single population",
    )
    demography_total.add_population(
        initial_size=rescale_population_size_by_selfing(
            population_size=params[0], selfing_rate=params[2], make_int=True
        ),
        description="Single population",
    )

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print(
            "".join(
                [
                    "\n",
                    "_" * 80,
                    "\nDemography for simulation\n",
                    "This demography is copy pasted from the creation of the",
                    "two phases, which we need to change the recombination rate",
                    "through time.\n\n",
                    str(demography_total),
                    "\n" + "=" * 80,
                ]
            ),
            file=logfile,
        )

    ts = simulate_transition_to_selfing(param_dict_list, rng)

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("simulated ancestry", file=logfile)

    seed = rng.integers(low=1, high=2**32 - 1, size=1)
    ts_mutated = msprime.sim_mutations(
        ts,
        rate=rule_parameters["mutrate"],
        random_seed=seed,
        model=None,
        start_time=None,
        end_time=None,
        discrete_genome=True,
        keep=True,
    )

    # log
    with open(log, "a", encoding="utf-8") as logfile:
        print(datetime.datetime.now(), end="\t", file=logfile)
        print("simulated mutations", file=logfile)

    sys.exit(
        "#" * 600
        + " inside simulate_treesequence_under_six_parameter_model\n"
        + " did you implement the six parameter model already?"
    )

    return ts_mutated


def simulate_transition_to_selfing(param_dict_list, rng):
    """Simulate transitions to selfing

    Technically we simulate here the six-parameter model with one change in
    population size and one change in selfing rate. This is a specified version
    of the piecewise constant populations size with piecewise constant
    recombination rate model.

    Args:
        param_dict_list: List with two dictionaries for the three phases of
            the simulation model
        rng: random number generator to produce the seeds for the simulations

    Returns:
        A single tree sequence
    """
    # ssimulate phase zero - selfing
    seed = rng.integers(low=1, high=2**32 - 1, size=1)
    ts_phase_zero = msprime.sim_ancestry(
        samples=param_dict_list[0]["samples"],
        demography=param_dict_list[0]["demography"],
        sequence_length=param_dict_list[0]["sequence_length"],
        discrete_genome=True,
        recombination_rate=param_dict_list[0]["recombination_rate_effective"],
        gene_conversion_rate=None,
        gene_conversion_tract_length=None,
        ploidy=2,  # The DTWF model only supports ploidy = 2
        model=param_dict_list[0]["msprime_model"],
        initial_state=None,
        start_time=None,
        end_time=param_dict_list[0]["end_time"],
        record_migrations=None,
        record_full_arg=None,
        num_labels=None,
        random_seed=seed,
        num_replicates=None,
        replicate_index=None,
        record_provenance=None,
    )

    # simulate phase one - outcrossing
    seed = rng.integers(low=1, high=2**32 - 1, size=1)
    ts_phase_one = msprime.sim_ancestry(
        demography=param_dict_list[1]["demography"],
        discrete_genome=True,
        recombination_rate=param_dict_list[1]["recombination_rate_effective"],
        gene_conversion_rate=None,
        gene_conversion_tract_length=None,
        ploidy=2,  # The DTWF model only supports ploidy = 2
        model=param_dict_list[1]["msprime_model"],
        initial_state=ts_phase_zero,
        start_time=param_dict_list[1]["start_time"],
        end_time=None,
        record_migrations=None,
        record_full_arg=None,
        num_labels=None,
        random_seed=seed,
        num_replicates=None,
        replicate_index=None,
        record_provenance=None,
    )

    # simplify tree sequence
    ts = ts_phase_one.simplify()

    return ts


def simulate_popsize_change_with_constant_selfing(param_dict, rng):
    """Simulate change in pop size

    Technically we simulate here the six-parameter model with one change in
    population size and one change in selfing rate. This is a specified version
    of the piecewise constant populations size with piecewise constant
    recombination rate model. This is the simulation of a piecewise constant
    population with constant selfing and only a single change in pop size.

    Args:
        param_dict_list: Dictionaries for the three phases of
            the simulation model
        rng: random number generator to produce the seeds for the simulations

    Returns:
        A single tree sequence
    """
    # ssimulate phase zero - selfing
    seed = rng.integers(low=1, high=2**32 - 1, size=1)
    ts_raw = msprime.sim_ancestry(
        samples=param_dict["samples"],
        demography=param_dict["demography"],
        sequence_length=param_dict["sequence_length"],
        discrete_genome=True,
        recombination_rate=param_dict["recombination_rate_effective"],
        gene_conversion_rate=None,
        gene_conversion_tract_length=None,
        ploidy=2,  # The DTWF model only supports ploidy = 2
        model=param_dict["msprime_model"],
        initial_state=None,
        start_time=None,
        end_time=None,
        record_migrations=None,
        record_full_arg=None,
        num_labels=None,
        random_seed=seed,
        num_replicates=None,
        replicate_index=None,
        record_provenance=None,
    )

    # simplify tree sequence
    ts = ts_raw.simplify()

    return ts


def simulate_transition_to_selfing_and_independent_change_of_pop_size(
    param_dict_list, rng
):
    sys.exit("#" * 600 + " inside simulate_treesequence_under_six_parameter_model")


def draw_parameter_from_prior(prior_definition, rng):
    """Draws parameters from a prior definition

    From a defined prior, we draw a single parameter value, this can be uniform,
    loguniform and returns float, or int. Ints will be drawn from a closed
    interval.

    Args:
        prior_definition: A list containing 4 values: min, max, numeric type
            (float, int), distribution type (uniform, loguniform)
        rng: The provided random number generator

    Returns:
        A single numeric value
    """
    assert prior_definition[2] in ["float", "int"], "unknown numeric format"
    assert prior_definition[3] in ["uniform", "loguniform"], "unknown distribution type"
    assert float(prior_definition[0]) <= float(
        prior_definition[1]
    ), "low value bigger than high value"

    # return the exact value if the boundaries are the same
    if float(prior_definition[0]) == float(prior_definition[1]):
        if prior_definition[2] == "int":
            drawn_parameter = int(float(prior_definition[0]))
        elif prior_definition[2] == "float":
            drawn_parameter = float(prior_definition[0])
        else:
            sys.exit(f"ERR: malformed prior definition")
    else:
        if prior_definition[2] == "float" and prior_definition[3] == "uniform":
            drawn_parameter = rng.uniform(
                low=float(prior_definition[0]),
                high=float(prior_definition[1]),
            )
        elif prior_definition[2] == "float" and prior_definition[3] == "loguniform":
            drawn_parameter = np.exp(
                rng.uniform(
                    low=np.log(float(prior_definition[0])),
                    high=np.log(float(prior_definition[1])),
                )
            )
        elif prior_definition[2] == "int" and prior_definition[3] == "uniform":
            drawn_parameter = rng.integers(
                low=int(float(prior_definition[0])),
                high=int(float(prior_definition[1])),
            )
        elif prior_definition[2] == "int" and prior_definition[3] == "loguniform":
            drawn_parameter = int(
                np.exp(
                    rng.uniform(
                        low=np.log(float(prior_definition[0])),
                        high=np.log(float(prior_definition[1])),
                    )
                )
            )
        else:
            sys.exit(f"ERR: malformed prior definition")

    return drawn_parameter


def rescale_recombination_rate_by_selfing(recombination_rate=1, selfing_rate=0):
    """Nordborgs (1997, 1999, 2000) coalescent with partial selfing

    Args:
        recombination_rate: float
        selfing_rate: self-fertilization rate

    Returns:
        float, rescaled recombination rate
    """
    assert selfing_rate >= 0, "selfing rate must be equal/bigger than zero"
    assert selfing_rate <= 1, "selfing rate must be smaller/equal than one"

    if selfing_rate == 1:
        warnings.warn(
            "You provided full selfing (100%); this results in zero recombination"
        )

    inbreeding_factor = selfing_rate / (2 - selfing_rate)
    rescaled_recombination_rate = recombination_rate * (1 - inbreeding_factor)
    return rescaled_recombination_rate


def rescale_population_size_by_selfing(
    population_size=1, selfing_rate=0, make_int=False
):
    """Nordborgs (1997, 1999, 2000) coalescent with partial selfing

    Args:
        population_size: float or int
        selfing_rate: self-fertilization rate
        make_int: boolean, whether to provide the rescaled population size as
            integer; default=False

    Returns:
        float, rescaled population size
    """
    assert selfing_rate >= 0, "selfing rate must be equal/bigger than zero"
    assert selfing_rate <= 1, "selfing rate must be smaller/equal than one"

    inbreeding_factor = selfing_rate / (2 - selfing_rate)
    rescaled_population_size = population_size * 1 / (1 + inbreeding_factor)

    if make_int:
        rescaled_population_size = int(rescaled_population_size)

    return rescaled_population_size


def find_breakpoints_for_TM_WIN(tsl, specs, rng, log):
    """Find good breakpoints for the discretization of TM_WIN diversity

    This function aims to equidistribute the information along the pairwise
    comparison of sequences into the number of discretized bins specified in
    specs. The discretization occurs data based and does not work on any
    theoretical expectations.

    Args:
        tsl: list of list of treesequences to take into account for the
            calculation of the breakpoints; e. g. if 100 reps and 5 loci:
            [[1, 2, ..., 5], ..., 100]
        specs: dict, specification how to calculate TM_WIN; expected to be a
            dictionary with at least following entries: winsize (window size of
            the sliding non-overlapping window), classes (number of discretized
            bins)

    Returns:
        np.array of floats, breakpoints to use
    """
    # get number of treeseqs
    num_trees = sum(1 for x in itertools.chain.from_iterable(tsl))

    # flatten
    tsl = itertools.chain.from_iterable(tsl)

    # loop through each ts and obtain the diversities
    collected_pairwise_diversities = []
    for tsid, ts in enumerate(tsl):
        # log
        with open(log, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                f"calculating diversity on {tsid+1}-th treeseq (of {num_trees})",
                file=logfile,
            )

        # chose one sample per individual
        sample_set = [
            this_individual.nodes[rng.integers(low=0, high=len(this_individual.nodes))]
            for this_individual in ts.individuals()
        ]

        # simplify the treesequence to the sampled haplotypes
        tshap = ts.simplify(samples=sample_set)
        assert (
            tshap.num_individuals == tshap.num_samples
        ), "you most likely did not sample one haplotype per individual"
        pairwise_diversities = tshap.diversity(
            sample_sets=list(itertools.combinations(tshap.samples(), 2)),
            windows=list(
                range(0, int(tshap.sequence_length), int(float(specs["winsize"])))
            )
            + [tshap.sequence_length],
            mode="site",
            span_normalise=False,
        )
        pairwise_diversities = list(itertools.chain.from_iterable(pairwise_diversities))
        collected_pairwise_diversities.append(pairwise_diversities)
    pairwise_diversities = np.array(
        list(itertools.chain.from_iterable(collected_pairwise_diversities))
    )
    del collected_pairwise_diversities

    # calculate which breakpoints to use
    breaks = [-0.5]  # first breakpoint must be left of first class (zero)
    unused_breaks = []
    freqs, counts = np.unique(pairwise_diversities, return_counts=True)
    for this_freq in range(int(pairwise_diversities.max())):
        proportion = (
            counts[(freqs > breaks[-1]) & (freqs <= this_freq)].sum() / counts.sum()
        )

        if proportion >= (1 / int(float(specs["classes"]))):
            breaks.append(this_freq + 0.5)
        else:
            unused_breaks.append(this_freq)

    # add more breakpoints if needed
    while len(breaks) < int(float(specs["classes"])):
        this_freq = unused_breaks[0]
        unused_breaks.remove(this_freq)
        breaks.append(this_freq + 0.5)
        if not len(unused_breaks):
            break

    breaks.append(np.Inf)
    breaks.sort()

    return np.array(breaks)


def calculate_ld_by_matrix(treeseq, npos, breaks, rng, log=False):
    """Calculate the ld on the treeseq

    This function aims to calculate a discretized LD pattern based on r2 for
    a set randomly chosen sites. This is very slow. However, if the treesequence
    provides the requesites to calculate the LD matrix, this is tried for at
    least the npos_threshold number of SNPs (default=1000). That is, if
    specs["npos"] is below, we try first to use the matrix approach with a
    subset of npos_threshold sites. If that fails, we use the snp times snp
    approach for the specified number of snps. If specs["npos"] is larger than
    the threshold, specs["npos"] will be used.

    Args:
        treeseq: treesequences
        specs: dict, specification how to calculate discretized LD; expected to
            be a dictionary with at least following entries: npos (number of
            SNPs to use in the calculation of an LD matrix)
            NOTE: breaks are not used
        breaks: breakpoints for discretization on the physical distance
        rng: random number generator to calculate which sites to use (see specs)
        log: str, logfile to print errors when calculating ld
        npos_threshold: int, number of sites to be used for the matrix LD calc
            (see tskit) approach (default=1000).

    Returns:
        np.array of floats being the discretized values for LD
    """
    try:
        # define site ids to calculate the ld on
        all_site_ids = list(site.id for site in treeseq.sites())
        if treeseq.num_sites > npos:
            siteids = rng.choice(all_site_ids, size=npos, replace=False, shuffle=False)
        else:
            siteids = all_site_ids

        # delete all sites that are not chosen
        site_ids_to_delete = [
            site_id for site_id in all_site_ids if not site_id in siteids
        ]

        treeseq = treeseq.delete_sites(site_ids_to_delete, record_provenance=True)

        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"deleted sites to calculate ld with matrix approach",
                    file=logfile,
                )

        # loop through all chose sites and calculate the ld_array on them
        ld_calculator = tskit.LdCalculator(treeseq)

        ld_matrix = ld_calculator.r2_matrix()

        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"calculated r2-matrix",
                    file=logfile,
                )

        sitepos = np.array([site.position for site in treeseq.sites()])
        ld_dist = np.array([abs(site - sitepos) for site in sitepos])

        r2_values = []
        r2_phys = []

        for r2_index, _ in np.ndenumerate(ld_matrix):
            if np.diff(r2_index) >= 0:
                r2_values.append(ld_matrix[r2_index])
                r2_phys.append(ld_dist[r2_index])

        r2_values = np.array(r2_values)
        r2_phys = np.array(r2_phys)

        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"calculated r2-distance-matrix",
                    file=logfile,
                )

        # discretize the r2 values by physical distance
        r2_classes = np.digitize(r2_phys, bins=breaks)
        r2_values = np.array(r2_values)

        final_r2_vector = []
        for classid in range(
            1, len(breaks) + 1
        ):  # because np.digitize results are 1-based
            if classid in r2_classes:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    final_r2_vector.append(np.nanmean(r2_values[r2_classes == classid]))
            else:
                final_r2_vector.append(0)

        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"discretized r2 values, returning final LD statistics",
                    file=logfile,
                )

        result = np.array(final_r2_vector)

    except:
        result = False

    return result


def calculate_ld(treeseq, specs, breaks, rng, log, npos_threshold=1000):
    """Calculate the ld on the treeseq

    This function aims to calculate a discretized LD pattern based on r2 for
    a set randomly chosen sites. This is very slow. However, if the treesequence
    provides the requesites to calculate the LD matrix, this is tried for at
    least the npos_threshold number of SNPs (default=1000). That is, if
    specs["npos"] is below, we try first to use the matrix approach with a
    subset of npos_threshold sites. If that fails, we use the snp times snp
    approach for the specified number of snps. If specs["npos"] is larger than
    the threshold, specs["npos"] will be used.

    Args:
        treeseq: treesequences
        specs: dict, specification how to calculate discretized LD; expected to
            be a dictionary with at least following entries: npos (number of
            SNPs to use in the calculation of an LD matrix)
            NOTE: breaks are not used
        breaks: breakpoints for discretization on the physical distance
        rng: random number generator to calculate which sites to use (see specs)
        log: str, logfile to print errors when calculating ld
        npos_threshold: int, number of sites to be used for the matrix LD calc
            (see tskit) approach (default=1000).

    Returns:
        np.array of floats being the discretized values for LD
    """
    npos = int(float(specs["npos"]))

    # calculate with the matrix approach, will return False if not successful
    ld_matrix_approach = calculate_ld_by_matrix(
        treeseq, max(npos, npos_threshold), breaks, rng, log=False
    )

    # check if matrix approach was successful, otherwise use the classical
    # but slow site x site approach
    if type(ld_matrix_approach) == bool and not ld_matrix_approach:
        # define site ids to calculate the ld on
        all_site_ids = list(site.id for site in treeseq.sites())
        if treeseq.num_sites > npos:
            siteids = rng.choice(all_site_ids, size=npos, replace=False, shuffle=False)
        else:
            siteids = all_site_ids

        # loop through all chose sites and calculate the ld_array on them
        ld_calculator = tskit.LdCalculator(treeseq)
        r2_values = []
        r2_phys = []
        ERRORS = False
        for site1, site2 in itertools.combinations(siteids, 2):
            try:
                r2_value = ld_calculator.r2(site1, site2)
                r2_values.append(r2_value)
                r2_phys.append(
                    abs(treeseq.site(site1).position - treeseq.site(site2).position)
                )
            except:
                ERRORS = True

        # log
        if ERRORS:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    "error with discrete genome calculating LD statistics",
                    file=logfile,
                )

        # discretize the r2 values by physical distance
        r2_classes = np.digitize(r2_phys, bins=breaks)
        r2_values = np.array(r2_values)

        final_r2_vector = []
        for classid in range(
            1, len(breaks) + 1
        ):  # because np.digitize results are 1-based
            if classid in r2_classes:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    final_r2_vector.append(np.nanmean(r2_values[r2_classes == classid]))
            else:
                final_r2_vector.append(0)
    else:
        final_r2_vector = ld_matrix_approach

    return np.array(final_r2_vector)


def calculate_tm_win(treeseq, specs, breaks):
    """Calculate the tm_win on the treeseq

    This function aims to calculate a discretized TM_WIN using the
    specifications provided

    Args:
        treeseq: treesequences; should be haploids
        specs: dict, specification how to calculate discretized LD; expected to
            be a dictionary with at least following entries: npos (number of
            SNPs to use in the calculation of an LD matrix), breaks (breakpoints
            how to discretize LD values based on physical distance of SNPs)

    Returns:
        np.array of floats being the discretized values for LD
    """
    # calculate pairwise diversities of all pairs
    pairwise_diversities = treeseq.diversity(
        sample_sets=list(itertools.combinations(treeseq.samples(), 2)),
        windows=list(
            range(0, int(treeseq.sequence_length), int(float(specs["winsize"])))
        )
        + [treeseq.sequence_length],
        mode="site",
        span_normalise=False,
    )

    # discretize the diversity (diversity into states)
    tm_win_classes = np.digitize(pairwise_diversities, bins=breaks) - 1

    # calculate the transition for each sequence of states
    tm_win_list = []
    for this_sequence in tm_win_classes.T:
        tm_win_list.append(
            transition_probability_matrix(this_sequence, len(breaks) - 1)
        )

    # calculate the mean for the transition probabilities
    tm_win_list = np.array(tm_win_list)
    tm_win_mean = tm_win_list.mean(axis=0)

    return np.array(list(itertools.chain.from_iterable(tm_win_mean)))


def transition_probability_matrix(sequence, nstates, logtransform=False):
    """From a sequence of states calculate the transition matrix

    The transition matrix provides the probability (frequentist) of states at
    positions n to n+1

    Args:
        sequence: np.array sequence of states
        nstates: integer, number of states

    Returns:
        np.matrix, transition probabilites between states
    """
    transition_matrix = np.zeros((nstates, nstates))

    for (i, j) in zip(sequence, sequence[1:]):
        transition_matrix[i][j] += 1

    # Convert to probabilities
    for row in transition_matrix:
        sum_row = sum(row)
        if sum_row > 0:
            if logtransform:
                row[:] = [np.log(sum_row) - np.log(element_row) for element_row in row]
            else:
                row[:] = [element_row / sum_row for element_row in row]

    return np.matrix(transition_matrix)


def filter_mask_for_disjoint_intervals(mask, log=False):
    """Filter mask for disjoint intervals

    tskit.delete_intervals() functions can take a mask argument. This function
    will provide a mask np.array(N, 2) that fits the requirements of the tskit
    function, i. e. all disjoint intervals

    Args:
        mask: np.array(N, 2), with N being the number of intervals
        log: path to log file or False

    Returns:
        np.array(N, 2)
    """
    # check for direction of interval, must be ascending
    for interval_id, (start, end) in enumerate(mask):
        if start > end:
            mask[interval_id] = np.array([end, start])

            # log
            if log:
                with open(log, "a", encoding="utf-8") as logfile:
                    print(datetime.datetime.now(), end="\t", file=logfile)
                    print(
                        f"changed interval direction for id (0-based): {interval_id}",
                        file=logfile,
                    )

    # check for disjointness
    overlapping_intervals = []
    for interval_id, (prev_end, start) in enumerate(zip(mask[:, 1], mask[1:, 0])):
        if prev_end > start:
            overlapping_intervals.append(interval_id)

    # fuse the overlapping intervals and repeat the loop
    counter = 0
    while len(overlapping_intervals):
        counter += 1
        if counter > 1e7:
            sys.exit(
                "#" * 600
                + " inside filter_mask_for_disjoint_intervals()\n"
                + "you passed the while-loop more than 1e7-times. Thus, I killed it.\n"
                + "In case you know that everything is correct, come here and change the \n"
                + "value to a higher threshold."
            )

        for overlapping_id in overlapping_intervals:
            mask[overlapping_id] = np.array(
                [mask[overlapping_id, 0], mask[overlapping_id + 1, 1]]
            )
            mask = np.delete(mask, overlapping_id + 1, axis=0)

            # re-index the overlapping id
            overlapping_intervals = [this_id - 1 for this_id in overlapping_intervals]

            # log
            if log:
                with open(log, "a", encoding="utf-8") as logfile:
                    print(datetime.datetime.now(), end="\t", file=logfile)
                    print(
                        f"fused intervals (0-based): {interval_id, interval_id+1}",
                        file=logfile,
                    )

        # create new list of intervals
        overlapping_intervals = []
        for interval_id, (prev_end, start) in enumerate(zip(mask[:, 1], mask[1:, 0])):
            if prev_end > start:
                overlapping_intervals.append(interval_id)

    # remove all intervals with zero length
    counter = 0
    while any(mask[:, 0] == mask[:, 1]):
        counter += 1
        if counter > 1e7:
            sys.exit(
                "#" * 600
                + " inside filter_mask_for_disjoint_intervals()\n"
                + "you passed the while-loop more than 1e7-times. Thus, I killed it.\n"
                + "In case you know that everything is correct, come here and change the \n"
                + "value to a higher threshold."
            )

        for maskid in np.where(mask[:, 0] >= mask[:, 1])[0]:
            mask = np.delete(mask, maskid, axis=0)

    # log
    if log:
        with open(log, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                "corrected the mask, zero-length intervals, filtered for disjoint intervals",
                file=logfile,
            )

    return mask


def check_masking_parameter(do_masking):
    """
    based on the parameter value of the rule, which supposedly is a boolean, but
    sometimes a string, return a boolean
    """
    #
    if type(do_masking) == str:
        do_masking = to_masking.casefold()
        assert do_masking in [
            "true",
            "false",
        ], "the masking parameter must be bool or str that is transferable to bool"

        # translate string to bool
        if do_masking == "true":
            do_masking = True
        elif do_maksing == "false":
            do_masking = False
        else:
            sys.exit(
                "#" * 600
                + "We are inside check_masking_parameter(), but under no condition we should ever reach here!"
            )

    else:
        assert (
            type(do_masking) == bool
        ), "the masking parameter must be bool or str that is transferable to bool"

    return do_masking


def random_subsets_from_iterable(
    iterable, size, nsam, rng, combination_threshold=1e5, log=False
):
    """Random subsamples from iterable

    This function randomly samples a subset of size size, repeats the sampling
    nsam times, but all being different. In precise, all possible combinations
    are indexed and from them will be sampled by random. If the number of
    possible combinations exceeds the combination_threshold, the subsamples
    will be random and possibly not unique. If the threshold is high, the
    calculation may take very long.

    Args:
        iterable: iterable to sample from
        size: subsample size
        nsam: number how often to repeat the sampling
        combination_threshold: numeric, threshold upon which the sampling will
            not take place indexbased, it may result in sample sets not being
            unique

    Returns:
        2d-np.arryay, sample sets
    """
    max_combinations = math.comb(len(iterable), min(size, len(iterable)))

    # limit sample size to max possible combinations
    nsam = min(nsam, max_combinations)

    # create sample subsets
    sampled_combinations = []

    # only define index to avoid same sample set if there are no huge number of simulations
    if max_combinations > 1e5:
        for sampelid in range(nsam):
            sampled_combinations.append(
                rng.choice(iterable, size=size, replace=False, shuffle=False)
            )

        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    "created random subsamples (no assurance for being unique)",
                    file=logfile,
                )

    else:
        indexes = rng.choice(
            range(max_combinations), size=nsam, replace=False, shuffle=False
        )
        max_idx = max(indexes)

        # generate samples
        sampled_combinations = []
        for combid, combo in zip(
            itertools.count(), itertools.combinations(iterable, size)
        ):
            if combid in indexes:
                sampled_combinations.append(combo)
            elif combid > max_idx:
                break

        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    "created index-based random (unique) subsamples",
                    file=logfile,
                )

    return np.array(sampled_combinations)


def create_subsets_from_treeseqlist(tsl, specs, rng, log=False):
    """Obtain subsets samples of treeseqs

    This function aims to provide subset samples of a treesequence list. The
    input is a treesequence list with one treesequence per region/chromosome
    with a large sample size. This function will resample from the treesequence
    and return a list of lists to maintain the region/chromosome structure of
    the input list.

    Args:
        tsl: treesequences
        specs: dict, two entries: 1) "num_observations": this is the
            maximum number of subsets to perform; 2) "nsam": sample size
        rng: random number generator to calculate which sites to use (see specs)
        log: str, logfile to print errors when calculating ld

    Returns:
        np.array of tskit.treesequences with added dimension compared the input
            tsl
    """
    # log
    if log:
        with open(log, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                "subsetting to final samplesize and resampling",
                file=logfile,
            )

    for treeid, (treeseq1, treeseq2) in enumerate(zip(tsl[:], tsl[1:])):
        assert list(treeseq1.samples()) == list(
            treeseq2.samples()
        ), f"samples different in different chromosomes: {treeid} and {treeid +1}"
        assert list(treeseq1.individuals()) == list(
            treeseq2.individuals()
        ), f"individuals different in different chromosomes: {treeid} and {treeid +1}"

    # log
    if log:
        with open(log, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                f"total samples/individuals: {tsl[0].num_samples}/{tsl[0].num_individuals}; subsample size {specs['nsam']}",
                file=logfile,
            )

    # assure the trees being haplotype
    tsl_haploid = []
    for tsid, ts in enumerate(tsl):
        # log
        if log:
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"haplodize {tsid+1}-th treeseq (of {len(tsl)})",
                    file=logfile,
                )

        # chose one sample per individual
        sample_set = [
            this_individual.nodes[rng.integers(low=0, high=len(this_individual.nodes))]
            for this_individual in ts.individuals()
        ]

        # simplify the treesequence to the sampled haplotypes
        tsl_haploid.append(ts.simplify(samples=sample_set))

    tsl = tsl_haploid
    del tsl_haploid

    # log
    if log:
        with open(log, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print(
                "trees are haploid",
                file=logfile,
            )

    # create sample list
    sample_set_list = random_subsets_from_iterable(
        tsl[0].samples(), specs["nsam"], specs["num_observations"], rng, log=log
    )

    # subset trees and provide multi-dim np.array of trees
    treeseq_list = [[] for _ in range(len(tsl))]
    for sample_id, sample_set in enumerate(sample_set_list, start=1):
        for tsid, ts in enumerate(tsl):
            treeseq_list[tsid].append(ts.simplify(samples=sample_set))

        # log
        if log and (not (sample_id % 50) or sample_id == len(sample_set_list)):
            with open(log, "a", encoding="utf-8") as logfile:
                print(datetime.datetime.now(), end="\t", file=logfile)
                print(
                    f"subsampled treeseqs {sample_id} of {len(sample_set_list)}",
                    file=logfile,
                )

    # log
    if log and (not (sample_id % 50) or sample_id == len(sample_set_list)):
        with open(log, "a", encoding="utf-8") as logfile:
            print(datetime.datetime.now(), end="\t", file=logfile)
            print("finished treesequence subsampling", file=logfile)

    return np.array(treeseq_list).T


def discretized_times(n=30, M=2):
    """Discretize time segments

    Here, we discretize the time segments following the quantiles of the exponential distribution. We determine the expected time to the first coalescence. The times are given in units 2*N_0, where N_0 is recommended to be fixed from Watterson's estimator.

    Args:
        n: int, number of segments
        M: int, number of haplotypes

    Returns:
        np.array with time segments
    """
    return np.array(
        [-np.log(1 - (i / n)) / (math.comb(M, 2)) for i in range(n)] + [np.Inf]
    )


def snp_freq_from_times(discretized_times, two_N_zero, mutrate, window_size):
    """Translate expected snps by provided times

    Calculate the expected number of SNPs per expected coalescent time. If boundaries are too interspersed, then we add the mid breaks between the smallest breaks allowing for entering another break.

    Args:
        discretized_times: np.array with discrete times, expect youngest time to be zero and largest to be np.Inf
        two_N_zero: float, the scaling parameter for the time, this value should usually be estimated on basis of a Watterson's estimator
        mutrate: float, per bp per gen mutation rate
        window_size: int (or float), to scale the window size
    """
    expected_snps = discretized_times * two_N_zero * window_size * mutrate

    # scale
    boundaries = [expected_snps[0]]
    for this_time in expected_snps[1:-1]:
        if int(this_time) == int(boundaries[-1]):
            continue
        else:
            boundaries.append(this_time)
    boundaries.append(expected_snps[-1])

    # delete too intersparsed
    delete_ids = []
    for bid, (boundary1, boundary2) in enumerate(
        zip(boundaries[:-1], boundaries[1:-1])
    ):
        if int(boundary1) + 1 == math.ceil(boundary2):
            delete_ids.append(bid + 1)
        else:
            continue

    filtered_boundaries = []
    for bid, boundary in enumerate(boundaries):
        if not bid in delete_ids:
            filtered_boundaries.append(boundary)
        else:
            continue
    boundaries = filtered_boundaries
    del filtered_boundaries
    boundaries.sort()

    # add more breaks if length is too small
    boundaries = list(boundaries)
    counter = 0
    while len(boundaries) < len(discretized_times):
        counter += 1
        for bid, (boundary1, boundary2) in enumerate(zip(boundaries, boundaries[1:])):
            if not np.isinf(boundary2):
                if math.ceil(boundary1) + 1 < int(boundary2):
                    boundaries.append(
                        math.ceil(boundary1)
                        + 0.5 * (math.ceil(boundary1) + int(boundary2))
                    )
                    boundaries.sort()
                    break
                else:
                    continue
            else:
                boundaries.append(boundary1 + 1)
                boundaries.sort()
        if counter < 1e8:
            continue
        else:
            assert False, "infinit loop"

    return np.array(boundaries)
