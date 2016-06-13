import Base.show
export ExpectationPropagation

type EPSite
    interface::Interface
    expectation_type::DataType
    expectation_entry::ScheduleEntry
    prior_schedule::Schedule
    prior::ProbabilityDistribution # payload of the incoming message, interface.partner.message always holds the cavity distribution

    function EPSite{T<:ProbabilityDistribution}(interface::Interface, expectation_type::Type{T})
        expectation_entry = convert(ScheduleEntry, interface, expectationRule!)
        expectation_entry.outbound_type = expectation_type
        new(interface, expectation_type, expectation_entry)
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
    execute::Function
    sites::Vector{EPSite}
    pre_schedule::Schedule
    post_schedule::Schedule
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
    #         execute site.prior_schedule
    #         expectation_dist = site.interface.message.payload
    #         if expectation_dist <: PartitionedDistribution:
    #             for every factor in expectation_dist:
    #                 calculate cavity distribution from prior + other factors
    #                 invoke expectationRule! on factor
    #         else:
    #             invoke expectationRule! (cavity distribution is equal to the result of algo.prior_schedule)
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
    site_interfaces = Set(keys(interface_to_site))
    dg = summaryDependencyGraph(graph)
    rdg = summaryDependencyGraph(graph, reverse_edges=true)
    influenced_by_site = [Set(children(site_interface, rdg)) for site_interface in site_interfaces]

    # Build prior schedules for all sites.
    # Dependencies that are independent of the sites end up in the pre-schedule,
    # which is only executed once at the start of the algorithm.
    pre_schedule = children(Interface[site_interface.partner for site_interface in site_interfaces], dg, breakers=site_interfaces)
    for i = 1:length(ep_sites)
        site = ep_sites[i]
        prev_site_idx = (i==1) ? length(ep_sites) : i-1
        message_types[site.interface] = site.expectation_type # Add sites to the fixed msg type dictionary
        # The prior schedule should include all messages that influence the prior on site i, and that are influenced by the previous site
        prior_schedule = children(site.interface.partner, dg, breakers=site_interfaces, restrict_to=influenced_by_site[prev_site_idx])
        ep_sites[i].prior_schedule = convert(Schedule, prior_schedule, sumProductRule!)
    end

    # build post-schedule
    post_schedule = children(targets, dg, breakers=union(site_interfaces, Set(pre_schedule)))

    # type inference
    pre_schedule = convert(Schedule, pre_schedule, sumProductRule!)
    post_schedule = convert(Schedule, post_schedule, sumProductRule!)
    inferMessageTypes!(pre_schedule, ep_sites, post_schedule, message_types)

    # # Build execute function
    # function exec(algorithm)
    #     # Init all sites with vague messages
    #     for site in algorithm.sites
    #         vague!(site.message.payload)
    #     end
    #     # Execute iterative schedule until stopping criterium is met
    #     for iteration_count = 1:algorithm.n_iterations
    #         execute(algorithm.iterative_schedule)
    #         # Check stopping criteria
    #         if algorithm.callback()
    #             break
    #         end
    #     end
    #     # Execute post convergence schedule once
    #     isempty(algorithm.post_convergence_schedule) || execute(algorithm.post_convergence_schedule)
    # end

    algo = ExpectationPropagation(graph,
                                  () -> true,
                                  ep_sites,
                                  pre_schedule,
                                  convert(Schedule, Interface[], sumProductRule!), # TODO: post-schedule
                                  n_iterations,
                                  callback)

    return algo
end

############################################
# Type inference and preparation
############################################

function inferMessageTypes!(pre_schedule::Schedule,
                            ep_sites::Vector{EPSite},
                            post_schedule::Schedule,
                            fixed_types::Dict{Interface,DataType})
    # Message type inference for all schedule entries in:
    #   - pre_schedule
    #   - site.prior_schedule for site in ep_sites
    #   - post_schedule
    # This fills the inbound_types and outbound_type fields of every ScheduleEntry

    schedule_entries = Dict{Interface, ScheduleEntry}() # lookup table
    fixed_type_interfaces = keys(fixed_types)
    #schedule_entries = vcat(pre_schedule, [site.prior_schedule for site in ep_sites]..., post_schedule)
end

function inferDistributionTypes!(   algo::ExpectationPropagation,
                                    recognition_distributions::Dict{Interface,DataType},
                                    message_types::Dict{Interface,DataType})
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}() # Lookup table from interface to schedule entry
    message_types_interfaces = keys(message_types)
    recognition_distributions_interfaces = keys(recognition_distributions)

    for entry in vcat(algo.iterative_schedule, algo.post_convergence_schedule)
        collectInboundTypes!(entry, schedule_entries, recognition_distributions, algo) # Fill entry.inbound_types
        outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
        if outbound_interface in message_types_interfaces
            setOutboundType!(entry, message_types[outbound_interface])
        elseif outbound_interface in recognition_distributions_interfaces
            setOutboundType!(entry, recognition_distributions[outbound_interface])
        end
        inferOutboundType!(entry) # Infer the outbound message type, or validate that there exists a suitable rule if the outbound type is already fixed
        schedule_entries[outbound_interface] = entry # Add entry to lookup table
    end

    return algo
end

function collectInboundTypes!{rule}(entry::ScheduleEntry{rule},
                                    schedule_entries::Dict{Interface, ScheduleEntry},
                                    recognition_distributions::Dict{Interface,DataType},
                                    algo::ExpectationPropagation)
    # Look up the types of the inbound messages for entry.
    # Fill entry.inbound_types
    entry.inbound_types = []

    for (id, interface) in enumerate(entry.node.interfaces)
        if (id == entry.outbound_interface_id) && (rule == SumProductRule)
            # Incoming msg on outbound interface is always Void for SumProductRule
            push!(entry.inbound_types, Void)
        elseif haskey(recognition_distributions, interface.partner)
            # Incoming msg from a site, so the type is given by the recognition distribution
            push!(entry.inbound_types, Message{recognition_distributions[interface.partner]})
        else
            # Incoming msg from earlier schedule entry
            push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
        end
    end

    return entry
end

function prepare!(algo::ExpectationPropagation)
    # Populate the graph with vague messages of the correct types
    for entry in vcat(algo.iterative_schedule, algo.post_convergence_schedule)
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedules
    compile!(algo.iterative_schedule, algo)
    compile!(algo.post_convergence_schedule, algo)

    return algo.graph.prepared_algorithm = algo
end

function compile!(entry::ScheduleEntry{ExpectationRule}, ::InferenceAlgorithm)
    # Generate entry.execute
    inbound_messages = [interface.partner.message for interface in entry.node.interfaces]

    return buildExecute!(entry, inbound_messages)
end
