import Base.show
export ExpectationPropagation

type EPSite
    interface::Interface
    expectation_type::DataType
    schedule::Schedule # last entry applies the ExpectationRule

    function EPSite{T<:ProbabilityDistribution}(interface::Interface, expectation_type::Type{T})
        new(interface, expectation_type)
    end
end

"""
Expectation propagation algorithm.

Usage:

    ExpectationPropagation(sites::Vector{Tuple{Interface, DataType}}; n_iterations, callback)
    ExpectationPropagation(target::Interface, sites::Vector{Tuple{Interface, DataType}}; n_iterations, callback)
    ExpectationPropagation(targets::Vector{Interface}, sites::Vector{Tuple{Interface, DataType}}; n_iterations, callback)
    ExpectationPropagation(graph::FactorGraph, sites::Vector{Tuple{Interface, DataType}}; n_iterations, callback)

Builds an expectation propagation message passing algorithm for the specified sites. Arguments:

    - `sites`: a list of `(interface, expectation_type)` tuples, where `interface` specifies where the expectation message
      is calculated, and `expectation_type <: ProbabilityDistribution` specifies the type of the expectation message.
    - `n_iterations`: the maximum number of iterations.
    - `callback`: a function that gets called after every iteration. If it returns true, the algorithm is terminated.
"""
type ExpectationPropagation <: InferenceAlgorithm
    graph::FactorGraph
    sites::Vector{EPSite}
    pre_schedule::Schedule                          # executed once before iterating of the sites
    post_schedule::Schedule                         # executed once after convergence
    n_iterations::Int64
    callback::Function
end

function show(io::IO, algo::ExpectationPropagation)
    println("ExpectationPropagation inference algorithm")
    println("    # sites: $(length(algo.sites))")
    println("    max. number of iterations: $(algo.n_iterations)")
    println("    callback function: $(algo.callback)")
end

function ExpectationPropagation(
            sites::Vector{Tuple{Interface, DataType}};
            kwargs...)

    ExpectationPropagation(Interface[], sites; kwargs...)
end

function ExpectationPropagation(
            target::Interface,
            sites::Vector{Tuple{Interface, DataType}};
            kwargs...)

    ExpectationPropagation([target], sites; kwargs...)
end

function ExpectationPropagation(
            graph::FactorGraph,
            sites::Vector{Tuple{Interface, DataType}};
            kwargs...)

    ExpectationPropagation(interfacesFacingWrapsOrBuffers(graph), sites; kwargs...)
end

function ExpectationPropagation(
            targets::Vector{Interface},
            sites::Vector{Tuple{Interface, DataType}};
            n_iterations::Int64 = 100,
            callback::Function = ( () -> false ),
            graph::FactorGraph=currentGraph(),
            message_types::Dict{Interface,DataType}=Dict{Interface,DataType}())
    # Algorithm pseudo code:
    #
    # init all sites with vague expectation messages.
    # execute algo.pre_schedule
    # while (algo.callback()!=true && algo.n_iterations not exceeded):
    #     for every site in algo.sites:
    #         execute site.schedule
    # execute algo.post_convergence_schedule

    (length(sites) > 0) || throw(ArgumentError("Specify at least one site"))

    # unzip sites
    interface_to_site = Dict{Interface, EPSite}()
    ep_sites = EPSite[]
    for (interface, expectation_type) in sites
        push!(ep_sites, EPSite(interface, expectation_type))
        interface_to_site[interface] = ep_sites[end]
    end

    # build schedules
    site_interfaces_set = Set(keys(interface_to_site))
    dg = summaryDependencyGraph(graph)
    rdg = summaryDependencyGraph(graph, reverse_edges=true)
    influenced_by_site = [Set(children(site_interface, rdg)) for site_interface in site_interfaces_set]

    # build pre-schedule
    site_dependencies = Interface[]
    for site in ep_sites
        push!(site_dependencies, site.interface.partner) # prior/cavity distribution
        for interface in neighbors(site.interface, dg)
            push!(site_dependencies, interface) # direct dependency of site
        end
    end
    pre_schedule = children(site_dependencies, dg, breakers=site_interfaces_set)

    # build site schedules
    # the last entry of site.schedule produces the expectation message
    for i = 1:length(ep_sites)
        site = ep_sites[i]
        prev_site_idx = (i==1) ? length(ep_sites) : i-1
        message_types[site.interface] = site.expectation_type # Add sites to the fixed msg type dictionary
        schedule = children([site.interface.partner; site.interface],
                            dg,
                            breakers = setdiff(site_interfaces_set, Set([site.interface])),
                            restrict_to = union(influenced_by_site[prev_site_idx], Set([site.interface])))
        site.schedule = ScheduleEntry[ScheduleEntry{SumProductRule}(interface) for interface in schedule[1:end-1]]
        expectation_entry = ScheduleEntry{ExpectationRule}(site.interface)
        expectation_entry.outbound_type = site.expectation_type
        push!(site.schedule, expectation_entry)
    end

    # build post-schedule
    post_schedule = children(targets, dg, breakers=union(site_interfaces_set, Set(pre_schedule)))

    # type inference
    pre_schedule = convert(Schedule, pre_schedule, SumProductRule)
    post_schedule = convert(Schedule, post_schedule, SumProductRule)
    for entry in pre_schedule
        inferTypes!(entry, message_types)
    end
    for site in ep_sites
        for entry in site.schedule
            inferTypes!(entry, message_types)
        end
    end
    for entry in post_schedule
        inferTypes!(entry, message_types)
    end

    algo = ExpectationPropagation(graph,
                                  ep_sites,
                                  pre_schedule,
                                  post_schedule,
                                  n_iterations,
                                  callback)

    return algo
end

############################################
# Type inference and preparation
############################################

function inferTypes!(   entry::ScheduleEntry{ExpectationRule},
                        inferred_outbound_types::Dict{Interface, DataType})
    # Collect inbound types
    entry.inbound_types = DataType[]
    for (id, interface) in enumerate(entry.node.interfaces)
        push!(entry.inbound_types, Message{inferred_outbound_types[interface.partner]})
    end

    # Infer outbound type
    outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
    if outbound_interface in keys(inferred_outbound_types)
        setOutboundType!(entry, inferred_outbound_types[outbound_interface])
    end
    inferOutboundType!(entry) # If entry.outbound_type is already set, this will validate that there is a suitable rule available
    inferred_outbound_types[outbound_interface] = entry.outbound_type

    return entry
end


function prepare!(algo::ExpectationPropagation)
    # Populate the graph with messages of the correct types and compile the schedules
    for entry in vcat(algo.pre_schedule, algo.post_schedule)
        interface = entry.node.interfaces[entry.outbound_interface_id]
        ensureMessage!(interface, entry.outbound_type)
    end
    for site in algo.sites
        for entry in site.schedule
            interface = entry.node.interfaces[entry.outbound_interface_id]
            ensureMessage!(interface, entry.outbound_type)
        end
    end

    # Compile the schedules
    compile!(algo.pre_schedule, algo)
    for site in algo.sites
        compile!(site.schedule, algo)
    end
    compile!(algo.post_schedule, algo)

    return algo.graph.prepared_algorithm = algo
end


function compile!(entry::ScheduleEntry{ExpectationRule}, ::InferenceAlgorithm)
    # Generate entry.execute
    inbound_messages = [interface.partner.message for interface in entry.node.interfaces]

    return buildExecute!(entry, inbound_messages)
end


function execute(algo::ExpectationPropagation)
    # Make sure algo is prepared
    (algo.graph.prepared_algorithm == algo) || prepare!(algo)

    # Init all sites with vague messages
    for site in algo.sites
        vague!(site.interface.message.payload)
    end

    # Execute pre schedule once
    isempty(algo.pre_schedule) || execute(algo.pre_schedule)

    # Loop over sites until stopping criterium is met
    for iteration_count = 1:algo.n_iterations
        for site in algo.sites
            execute(site.schedule)
        end

        # Check stopping criteria
        if algo.callback()
            break
        end
    end

    # Execute post convergence schedule once
    isempty(algo.post_schedule) || execute(algo.post_schedule)
end