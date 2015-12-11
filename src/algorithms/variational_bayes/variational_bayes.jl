import Base.show
export VariationalBayes

include("subgraph.jl")
include("q_distribution.jl")
include("q_factorization.jl")
include("scheduler.jl")
include("message_passing.jl")

type VariationalBayes <: InferenceAlgorithm
    execute::Function
    factorization::QFactorization
    q_distributions::Dict{Tuple{Node,Subgraph},QDistribution}
    n_iterations::Int64
end

function VariationalBayes(graph::FactorGraph=currentGraph(); n_iterations::Int64=50)
    # Generates a VariationalBayes algorithm that propagates messages to all write buffers and wraps.
    # Uses a mean field factorization and autoscheduler
    factorization = factorize(graph) # Mean field factorization
    generateVariationalBayesSchedule!(factorization, graph) # Generate and store internal and external schedules on factorization subgraphs
    q_distributions = vagueQDistributions(factorization) # Initialize vague q distributions

    function exec(algorithm)
        vagueQDistributions!(algorithm.q_distributions) # Reset q distributions before next step
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.factorization, algorithm.q_distributions) # For all subgraphs, execute internal and external schedules
        end
    end

    algo = VariationalBayes(exec, factorization, q_distributions, n_iterations)

    for factor in algo.factorization.factors
        compile!(factor.subgraph.internal_schedule, algo)
    end

    return algo
end

function VariationalBayes(graph::FactorGraph, cluster_edges...; n_iterations::Int64=50)
    # Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on graph
    # Uses a structured factorization with elements in cluster_edges defining separate subgraphs
    factorization = QFactorization(graph)
    for cluster in cluster_edges
        factorization = factorize!(cluster, factorization) # Structured factorization with elements in cluster_edges defining separate subgraphs
    end
    generateVariationalBayesSchedule!(factorization, graph) # Generate and store internal and external schedules on factorization subgraphs
    q_distributions = vagueQDistributions(factorization) # Initialize vague q distributions

    function exec(algorithm)
        vagueQDistributions!(algorithm.q_distributions) # Reset q distributions before next step
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.factorization, algorithm.q_distributions) # For all subgraphs, execute internal and external schedules
        end
    end

    algo = VariationalBayes(exec, factorization, q_distributions, n_iterations)
    for factor in algo.factorization.factors
        compile!(factor.subgraph.internal_schedule, algo)
    end

    return algo
end
VariationalBayes(cluster_edges...; n_iterations::Int64=50) = VariationalBayes(currentGraph(), cluster_edges...; n_iterations=n_iterations)

function compile!(schedule_entry::ScheduleEntry, ::Type{Val{symbol("ForneyLab.vmp!")}}, algo::VariationalBayes)
    # Compile ScheduleEntry objects for SumProduct algorithm
    # Generates schedule_entry.execute function

    # Collect references to all required inbound messages for executing message computation rule
    outbound_interface = schedule_entry.interface
    outbound_interface_index = 0
    inbounds = Any[]
    for j = 1:length(outbound_interface.node.interfaces)
        interface = outbound_interface.node.interfaces[j]
        if is(interface, outbound_interface)
            # We don't need the inbound message on the outbound interface
            outbound_interface_index = j
            push!(inbounds, nothing) # This interface is outbound, push "nothing"
        else
            if !haskey(algo.factorization.edge_to_subgraph, interface.edge) || !haskey(algo.factorization.edge_to_subgraph, outbound_interface.edge)
                # Inbound and/or outbound edge is not explicitly listed in the algo.fields.
                # This is possible if one of those edges is internal to a composite node.
                # We will default to sum-product message passing, and consume the message on the inbound interface.
                # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph, so it is safe to assume sum-product.
                try push!(inbounds, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
                break
            end

            # Should we require the inbound message or marginal?
            if is(algo.factorization.edge_to_subgraph[interface.edge], algo.factorization.edge_to_subgraph[outbound_interface.edge])
                # Both edges in same subgraph, require message
                try push!(inbounds, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
            else
                # A subgraph border is crossed, require marginal
                # The factor is the set of internal edges that are in the same subgraph
                sg = algo.factorization.edge_to_subgraph[interface.edge]
                push!(inbounds, algo.q_distributions[(outbound_interface.node, sg)].distribution)
            end
        end
    end

    # Assign the "compiled" update rule as an anomynous function to the schedule entry execute field
    schedule_entry.execute = ( () -> vmp!(node, outbound_interface_index, inbounds...) )

    return schedule_entry
end
