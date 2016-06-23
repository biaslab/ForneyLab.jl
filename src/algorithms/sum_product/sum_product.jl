import Base.show
export SumProduct

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
                    kwargs...)
    return SumProduct(interfacesFacingWrapsOrBuffers(graph); graph=graph, kwargs...)
end

function SumProduct(outbound_interface::Interface;
                    kwargs...)
    return SumProduct(Interface[outbound_interface]; kwargs...)
end

function SumProduct(edge::Edge;
                    message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                    graph::FactorGraph=currentGraph())
    # Generate schedule
    dg = summaryDependencyGraph(graph)
    schedule = convert(Schedule, children([edge.head, edge.tail], dg), SumProductRule)

    # Infer message types
    for entry in schedule
        inferTypes!(entry, message_types)
    end

    function exec(algorithm)
        execute(algorithm.schedule)
        calculateMarginal!(edge)
    end

    return SumProduct(graph, exec, schedule)
end

function SumProduct(outbound_interfaces::Vector{Interface};
                    message_types::Dict{Interface,DataType}=Dict{Interface,DataType}(),
                    graph::FactorGraph=currentGraph())
    # Generate schedule
    dg = summaryDependencyGraph(graph)
    schedule = convert(Schedule, children(outbound_interfaces, dg), SumProductRule)

    # Infer message types
    for entry in schedule
        inferTypes!(entry, message_types)
    end

    exec(algorithm) = execute(algorithm.schedule)

    return SumProduct(graph, exec, schedule)
end

############################################
# Type inference and preparation
############################################

function inferTypes!(   entry::ScheduleEntry{SumProductRule},
                        inferred_outbound_types::Dict{Interface, DataType})
    # Collect inbound types
    entry.inbound_types = DataType[]
    for (id, interface) in enumerate(entry.node.interfaces)
        if id == entry.outbound_interface_id
            push!(entry.inbound_types, Void)
        else
            push!(entry.inbound_types, Message{inferred_outbound_types[interface.partner]})
        end
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

function prepare!(algo::SumProduct)
    # Populate the graph with vague messages of the correct types
    for entry in algo.schedule
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedule (define entry.execute)
    compile!(algo.schedule, algo)

    return algo.graph.prepared_algorithm = algo
end

function compile!(entry::ScheduleEntry{SumProductRule}, ::InferenceAlgorithm)
    # Generate entry.execute

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
