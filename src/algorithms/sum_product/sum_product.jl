import Base.show
export SumProduct

include("scheduler.jl")

abstract AbstractSumProduct <: InferenceAlgorithm

"""
Sum-product message passing algorithm.

Usage:

    SumProduct(graph::Graph; post_processing_functions)
    SumProduct(outbound_interface::Interface; post_processing_functions)
    SumProduct(partial_list::Vector{Interface}; post_processing_functions)
    SumProduct(edge::Edge; post_processing_functions)
"""
type SumProduct <: AbstractSumProduct
    graph::FactorGraph
    execute::Function
    schedule::Schedule
end

function show(algo::SumProduct)
    println("SumProduct inference algorithm")
    println("    message passing schedule length: $(length(algo.schedule))")
    println("Use show(algo.schedule) to view the message passing schedule.")
end


############################################
# SumProduct algorithm constructors
############################################

function SumProduct(graph::FactorGraph=currentGraph(); post_processing_functions=Dict{Interface, Function}())
    schedule = generateSumProductSchedule(graph)
    setPostProcessing!(schedule, post_processing_functions)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function SumProduct(outbound_interface::Interface; post_processing_functions=Dict{Interface, Function}(), graph::FactorGraph=currentGraph())
    schedule = generateSumProductSchedule(outbound_interface)
    setPostProcessing!(schedule, post_processing_functions)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function SumProduct(partial_list::Vector{Interface}; post_processing_functions=Dict{Interface, Function}(), graph::FactorGraph=currentGraph())
    schedule = generateSumProductSchedule(partial_list)
    setPostProcessing!(schedule, post_processing_functions)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function SumProduct(edge::Edge; post_processing_functions=Dict{Interface, Function}(), graph::FactorGraph=currentGraph())
    schedule = generateSumProductSchedule([edge.head, edge.tail])
    setPostProcessing!(schedule, post_processing_functions)
    function exec(algorithm)
        execute(algorithm.schedule)
        calculateMarginal!(edge)
    end

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end


############################################
# Type inference and preparation
############################################

function inferDistributionTypes!(algo::AbstractSumProduct)
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}()

    for entry in algo.schedule
        collectInboundTypes!(entry, schedule_entries, algo) # SumProduct specific collection of inbound types
        inferOutboundType!(entry) # For the SumProduct algorithm, the only allowed update rule is the sumProductRule! rule

        outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
        schedule_entries[outbound_interface] = entry # Assign schedule entry to lookup dictionary
    end

    return algo
end

function collectInboundTypes!(entry::ScheduleEntry, schedule_entries::Dict{Interface, ScheduleEntry}, ::SumProduct)
    # Infers inbound types for node relative to the outbound interface
    entry.inbound_types = []

    for (id, interface) in enumerate(entry.node.interfaces)
        if id == entry.outbound_interface_id
            push!(entry.inbound_types, Void)
        else
            push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
        end
    end

    return entry
end

function prepare!(algo::SumProduct)
    # Populate the graph with vague messages of the correct types
    for entry in algo.schedule
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedule (define entry.execute)
    compile!(algo.schedule, algo)

    return algo.graph.prepared_algorithm = algo
end

function compile!(entry::ScheduleEntry, ::Type{Val{symbol(sumProductRule!)}}, ::InferenceAlgorithm)
    # Generate entry.execute for schedule entry with sumProductRule! update rule

    inbound_rule_arguments = []
    # Add inbound messages to inbound_rule_arguments
    for (id, interface) in enumerate(entry.node.interfaces)
        if id == entry.outbound_interface_id
            push!(inbound_rule_arguments, nothing)
        else
            push!(inbound_rule_arguments, interface.partner.message)
        end
    end

    return buildExecute!(entry, inbound_rule_arguments)
end
