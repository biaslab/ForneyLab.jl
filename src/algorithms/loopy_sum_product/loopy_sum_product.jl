import Base.show
export LoopySumProduct

"""
Loopy sum-product message passing algorithm.

You should provide breaker messages in every cycle in the graph.
If you don't, a cycle error will be thrown.

Usage:

    LoopySumProduct(graph::Graph; breaker_messages, message_types, n_iterations)
    LoopySumProduct(outbound_interface::Interface; breaker_messages, message_types, n_iterations)
    LoopySumProduct(outbound_interfaces::Vector{Interface}; breaker_messages, message_types, n_iterations)
"""
type LoopySumProduct <: AbstractSumProduct
    graph::FactorGraph
    pre_schedule::Schedule
    iterative_schedule::Schedule
    post_schedule::Schedule
    breaker_messages::Dict{Interface, Message} # Sites for breaker message initializations
    n_iterations::Int64
end

function show(io::IO, algo::LoopySumProduct)
    println("LoopySumProduct inference algorithm")
    println("    pre-schedule length: $(length(algo.pre_schedule))")
    println("    iterative schedule length: $(length(algo.iterative_schedule))")
    println("    post-schedule length: $(length(algo.post_schedule))")
    println("    number of iterations: $(algo.n_iterations)")
    println("Use show(algo.pre_schedule), show(algo.iterative_schedule), show(algo.post_schedule) to view the message passing schedules.")
end

############################################
# LoopySumProduct algorithm constructors
############################################

function LoopySumProduct(   graph::FactorGraph=currentGraph();
                            args...)
    return LoopySumProduct(interfacesFacingWrapsOrBuffers(graph); graph=graph, args...)
end

function LoopySumProduct(   outbound_interface::Interface;
                            args...)
    return LoopySumProduct(Interface[outbound_interface]; args...)
end

function LoopySumProduct(   outbound_interfaces::Vector{Interface};
                            breaker_messages=Dict{Interface, Message}(),
                            message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                            n_iterations=50,
                            graph::FactorGraph=currentGraph())
    !isempty(breaker_messages) || error("You should specify at least one breaker message")

    # Generate schedules
    pre_schedule = Interface[]          # all messages that are required for message passing in cycles, but that are not part of a cycle themselves
    iterative_schedule = Interface[]    # all required messages that are part of a cycle
    post_schedule = Interface[]         # all required messages that are not part of pre_schedule and iterative_schedule
    dg = summaryDependencyGraph(graph)  # dependency graph used for scheduling
    rdg = summaryDependencyGraph(graph, reverse_edges=true) # reverse dependency graph used for checking which messages depend on a given message
    breaker_sites = collect(keys(breaker_messages))
    # build pre_schedule and iterative_schedule
    influenced_by_breakers = children(breaker_sites, rdg, allow_cycles=true)
    for interface in children(breaker_sites, dg, allow_cycles=true)
        if interface in influenced_by_breakers
            push!(iterative_schedule, interface)
        else
            push!(pre_schedule, interface)
        end
    end
    # build post_schedule
    for interface in children(outbound_interfaces, dg, breakers=union(Set(pre_schedule), Set(iterative_schedule)))
        push!(post_schedule, interface)
    end
    pre_schedule = convert(Schedule, pre_schedule, SumProductRule)
    iterative_schedule = convert(Schedule, iterative_schedule, SumProductRule)
    post_schedule = convert(Schedule, post_schedule, SumProductRule)

    # Infer message types
    for (interface, message) in breaker_messages
        message_types[interface] = typeof(message.payload)
    end
    # pre_schedule type inference
    for entry in pre_schedule
        inferTypes!(entry, message_types)
    end
    # iterative_schedule type inference
    breaker_entries = ScheduleEntry[]
    for entry in iterative_schedule
        interface = entry.node.interfaces[entry.outbound_interface_id]
        if interface in breaker_sites
            push!(breaker_entries, entry) # inbounds for the breaker entry are not known yet
            continue
        else
            inferTypes!(entry, message_types)
        end
    end
    for entry in breaker_entries
        inferTypes!(entry, message_types)
    end
    # post_schedule type inference
    for entry in post_schedule
        inferTypes!(entry, message_types)
    end

    return LoopySumProduct( graph,
                            pre_schedule,
                            iterative_schedule,
                            post_schedule,
                            breaker_messages,
                            n_iterations)
end


############################################
# Type inference and breaker message reset
############################################

function prepare!(algo::LoopySumProduct)
    # Populate the graph with vague messages of the correct types
    for schedule in (algo.pre_schedule, algo.iterative_schedule, algo.post_schedule)
        for entry in schedule
            ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
        end
    end

    # Compile the schedules (define entry.execute)
    compile!(algo.pre_schedule, algo)
    compile!(algo.iterative_schedule, algo)
    compile!(algo.post_schedule, algo)

    return algo.graph.prepared_algorithm = algo
end

function resetBreakerMessages(algo::LoopySumProduct)
    # Before starting a new run, the breaker messages should be reset to their initial values
    for (interface, source_message) in algo.breaker_messages
        # Populate the original distribution parameters on interface with the breaker parameters
        injectParameters!(interface.message.payload, source_message.payload)
    end

    return algo
end

function execute(algo::LoopySumProduct)
    (algo.graph.prepared_algorithm == algo) || prepare!(algo)

    resetBreakerMessages(algo)

    isempty(algo.pre_schedule) || execute(algo.pre_schedule)
    for iteration = 1:algo.n_iterations
        execute(algo.iterative_schedule)
    end
    isempty(algo.post_schedule) || execute(algo.post_schedule)
end
