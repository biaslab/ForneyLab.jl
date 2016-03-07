import Base.show
export LoopySumProduct

"""
Loopy sum-product message passing algorithm.

Usage:

    LoopySumProduct(graph::Graph; breaker_messages, message_types, n_iterations)
    LoopySumProduct(outbound_interface::Interface; breaker_messages, message_types, n_iterations)
"""
type LoopySumProduct <: AbstractSumProduct
    graph::FactorGraph
    execute::Function
    schedule::Schedule
    breaker_messages::Dict{Interface, Message} # Sites for breaker message initializations
    n_iterations::Int64
end

function show(algo::LoopySumProduct)
    println("LoopySumProduct inference algorithm")
    println("    message passing schedule length: $(length(algo.schedule))")
    println("    number of iterations: $(algo.n_iterations)")
    println("Use show(algo.schedule) to view the message passing schedule.")
end

############################################
# LoopySumProduct algorithm constructors
############################################

function LoopySumProduct(   graph::FactorGraph=currentGraph();
                            breaker_messages=Dict{Interface, Message}(),
                            message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                            n_iterations=50)
    # Generates a LoopySumProduct algorithm that propagates messages to all wraps and write buffers.

    schedule = generateSumProductSchedule(graph, breaker_sites=Set(keys(breaker_messages)))

    function exec(algorithm)
        resetBreakerMessages(algorithm)
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.schedule)
        end
    end

    algo = LoopySumProduct(graph, exec, schedule, breaker_messages, n_iterations)
    inferDistributionTypes!(algo, message_types)

    return algo
end

function LoopySumProduct(   outbound_interface::Interface;
                            breaker_messages=Dict{Interface, Message}(),
                            message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                            n_iterations=50,
                            graph::FactorGraph=currentGraph())
    # Generates a LoopySumProduct algorithm to calculate the outbound message on outbound_interface.

    schedule = generateSumProductSchedule(outbound_interface, breaker_sites=Set(keys(breaker_messages)))

    function exec(algorithm)
        resetBreakerMessages(algorithm)
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.schedule)
        end
    end

    algo = LoopySumProduct(graph, exec, schedule, breaker_messages, n_iterations)
    inferDistributionTypes!(algo, message_types)

    return algo
end


############################################
# Type inference and breaker message reset
############################################

# Most methods are inherited from SumProduct

function collectInboundTypes!(entry::ScheduleEntry, schedule_entries::Dict{Interface, ScheduleEntry}, algo::LoopySumProduct)
    # Infers inbound types for node relative to the outbound interface
    entry.inbound_types = []

    for (id, interface) in enumerate(entry.node.interfaces)
        if id == entry.outbound_interface_id
            push!(entry.inbound_types, Void)
        elseif haskey(algo.breaker_messages, interface.partner) # A breaker message is pre-set on the partner interface, push breaker message type
            push!(entry.inbound_types, typeof(algo.breaker_messages[interface.partner]))
        else
            push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
        end
    end

    return entry
end

function prepare!(algo::LoopySumProduct)
    # Populate the graph with vague messages of the correct types
    for (interface, message) in algo.breaker_messages # Preset breaker messages
        ensureMessage!(interface, typeof(message.payload))
    end

    for entry in algo.schedule # Preset other messages
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedule (define entry.execute)
    compile!(algo.schedule, algo)

    return algo.graph.prepared_algorithm = algo
end

function resetBreakerMessages(algo::LoopySumProduct)
    # Before starting a new iteration, the breaker messages should be reset to their initial values
    for (interface, source_message) in algo.breaker_messages
        # Populate the original distribution parameters on interface with the breaker parameters
        injectParameters!(interface.message.payload, source_message.payload)
    end

    return algo
end
