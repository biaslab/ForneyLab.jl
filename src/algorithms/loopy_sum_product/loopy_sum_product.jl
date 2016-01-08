export LoopySumProduct

type LoopySumProduct <: AbstractSumProduct
    execute::Function
    schedule::Schedule
    breaker_messages::Dict{Interface, Message} # Sites for breaker message initializations
    n_iterations::Int64
end


############################################
# LoopySumProduct algorithm constructors
############################################

function LoopySumProduct(graph::FactorGraph=currentGraph(); post_processing_functions=Dict{Interface, Function}(), breaker_messages=Dict{Interface, Message}(), n_iterations=50)
    # Generates a LoopySumProduct algorithm that propagates messages to all wraps and write buffers.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(graph, breaker_sites=Set(keys(breaker_messages)))
    setPostProcessing!(schedule, post_processing_functions)

    function exec(algorithm)
        resetBreakerMessages(algorithm)
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.schedule)
        end
    end

    algo = LoopySumProduct(exec, schedule, breaker_messages, n_iterations)
    inferDistributionTypes!(algo)

    return algo
end

function LoopySumProduct(outbound_interface::Interface; post_processing_functions=Dict{Interface, Function}(), breaker_messages=Dict{Interface, Message}(), n_iterations=50)
    # Generates a LoopySumProduct algorithm to calculate the outbound message on outbound_interface.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(outbound_interface, breaker_sites=Set(keys(breaker_messages)))
    setPostProcessing!(schedule, post_processing_functions)

    function exec(algorithm)
        resetBreakerMessages(algorithm)
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.schedule)
        end
    end

    algo = LoopySumProduct(exec, schedule, breaker_messages, n_iterations)
    inferDistributionTypes!(algo)

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
            push!(entry.inbound_types, Void) # Outbound interface, push Void
        else
            if haskey(algo.breaker_messages, interface.partner) # A breaker message is pre-set on the partner interface, push message type
                push!(entry.inbound_types, typeof(algo.breaker_messages[interface.partner]))
            else
                push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
            end
        end
    end

    return entry
end

function resetBreakerMessages(algo::LoopySumProduct)
    # Before starting a new iteration, the breaker messages should be reset to their initial values
    for (interface, source_message) in algo.breaker_messages
        # Populate the original distribution parameters on interface with the breaker parameters
        injectParameters!(interface.message.payload, source_message.payload)
    end

    return algo
end

