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

#------------------------------
# VariationalBayes constructors
#------------------------------

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

    global current_algorithm = VariationalBayes(exec, factorization, q_distributions, n_iterations)
    return current_algorithm
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

    global current_algorithm = VariationalBayes(exec, factorization, q_distributions, n_iterations)
    return current_algorithm
end
VariationalBayes(cluster_edges...; n_iterations::Int64=50) = VariationalBayes(currentGraph(), cluster_edges...; n_iterations=n_iterations)


#---------------------------------------------------
# Construct algorithm specific update-call signature
#---------------------------------------------------

function collectInbounds(outbound_interface::Interface, ::Type{Val{symbol("ForneyLab.vmp!")}})
    # VMP specific method to collect all required inbound messages and marginals in an array.
    # This array is used to call the node update function (vmp!)
    # outbound_interface: the interface on which the outbound message will be updated
    # Returns: (outbound interface id, array of inbound messages and marginals)

    outbound_interface_index = 0
    inbounds = Array(Any, 0)
    for j = 1:length(outbound_interface.node.interfaces)
        interface = outbound_interface.node.interfaces[j]
        if is(interface, outbound_interface)
            # We don't need the inbound message on the outbound interface
            outbound_interface_index = j
            push!(inbounds, nothing) # This interface is outbound, push "nothing"
        else
            if !haskey(ForneyLab.current_algorithm.factorization.edge_to_subgraph, interface.edge) || !haskey(ForneyLab.current_algorithm.factorization.edge_to_subgraph, outbound_interface.edge)
                # Inbound and/or outbound edge is not explicitly listed in the ForneyLab.current_algorithm.fields.
                # This is possible if one of those edges is internal to a composite node.
                # We will default to sum-product message passing, and consume the message on the inbound interface.
                # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph, so it is safe to assume sum-product.
                try return push!(inbounds, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
            end

            # Should we require the inbound message or marginal?
            if is(ForneyLab.current_algorithm.factorization.edge_to_subgraph[interface.edge], ForneyLab.current_algorithm.factorization.edge_to_subgraph[outbound_interface.edge])
                # Both edges in same subgraph, require message
                try push!(inbounds, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
            else
                # A subgraph border is crossed, require marginal
                # The factor is the set of internal edges that are in the same subgraph
                sg = ForneyLab.current_algorithm.factorization.edge_to_subgraph[interface.edge]
                push!(inbounds, ForneyLab.current_algorithm.q_distributions[(outbound_interface.node, sg)].distribution)
            end
        end
    end

    return (outbound_interface_index, inbounds)
end