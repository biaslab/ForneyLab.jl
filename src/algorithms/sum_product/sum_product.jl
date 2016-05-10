import Base.show
export SumProduct

include("scheduler.jl")

abstract AbstractSumProduct <: InferenceAlgorithm

"""
Sum-product message passing algorithm.

Usage:

    SumProduct(graph::Graph)
    SumProduct(outbound_interface::Interface)
    SumProduct(partial_list::Vector{Interface})
    SumProduct(edge::Edge)

Optionally, keyword argument `message_types::Dict{Interface,Any}` can used to constain messages to a specific distribution type.
To force the use of a specific approximation method, use a tuple: `(dist_type::DataType, approximation::Symbol)`. Examples:

    message_types = Dict{Interface,DataType}(my_node.i[:out] => Gaussian)
    message_types = Dict{Interface,DataType}(my_node.i[:out] => Approximation{Gaussian,:laplace})
"""
type SumProduct <: AbstractSumProduct
    graph::FactorGraph
    execute::Function
    schedule::Schedule
end

function show(io::IO, algo::SumProduct)
    println("SumProduct inference algorithm")
    println("    message passing schedule length: $(length(algo.schedule))")
    println("Use show(algo.schedule) to view the message passing schedule.")
end


############################################
# SumProduct algorithm constructors
############################################

function SumProduct(graph::FactorGraph=currentGraph();
                    message_types::Dict{Interface,DataType}=Dict{Interface,DataType}())
    schedule = generateSumProductSchedule(graph)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo, message_types)

    return algo
end

function SumProduct(outbound_interface::Interface;
                    message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                    graph::FactorGraph=currentGraph())
    schedule = generateSumProductSchedule(outbound_interface)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo, message_types)

    return algo
end

function SumProduct(partial_list::Vector{Interface};
                    message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                    graph::FactorGraph=currentGraph())
    schedule = generateSumProductSchedule(partial_list)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo, message_types)

    return algo
end

function SumProduct(edge::Edge;
                    message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                    graph::FactorGraph=currentGraph())
    schedule = generateSumProductSchedule([edge.head, edge.tail])
    function exec(algorithm)
        execute(algorithm.schedule)
        calculateMarginal!(edge)
    end

    algo = SumProduct(graph, exec, schedule)
    inferDistributionTypes!(algo, message_types)

    return algo
end


############################################
# Type inference and preparation
############################################

function inferDistributionTypes!(algo::AbstractSumProduct, message_types::Dict{Interface,DataType})
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}()

    for entry in algo.schedule
        collectInboundTypes!(entry, schedule_entries, algo) # SumProduct specific collection of inbound types
        outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
        if outbound_interface in keys(message_types)
            setOutboundType!(entry, message_types[outbound_interface])
        end
        inferOutboundType!(entry) # If the outbound type is fixed, this will check if there is a suitable rule available

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
