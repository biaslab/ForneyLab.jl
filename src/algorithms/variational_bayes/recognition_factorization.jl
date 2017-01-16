import Base.factor
export RecognitionFactorization, currentRecognitionFactorization
export show, factor, factorizeMeanField, initialize

type RecognitionFactorization
    graph::FactorGraph

    # Recognition factorization
    subgraphs::Dict{Symbol, Subgraph}
    subgraph_id_counter::Int64 # Counter for next default subgraph identifier
    edge_to_subgraph::Dict{Edge, Subgraph}
    node_subgraph_to_internal_edges::Dict{Tuple{Node, Subgraph}, Set{Edge}}

    # Recognition distributions
    initial_recognition_distributions::Partitioned
    recognition_distributions::Partitioned
    node_subgraph_to_recognition_distribution::Dict{Tuple{Node, Subgraph}, ProbabilityDistribution}

    # This update schedule gets copied to the algorithm as the default
    default_update_order::Vector{Subgraph}

    RecognitionFactorization(graph) = new(  graph,
                                            Dict{Symbol, Subgraph}(),
                                            0,
                                            Dict{Edge, Subgraph}(),
                                            Dict{Tuple{Node, Subgraph}, Set{Edge}}(),
                                            Partitioned{ProbabilityDistribution,0}(ProbabilityDistribution[]),
                                            Partitioned{ProbabilityDistribution,0}(ProbabilityDistribution[]),
                                            Dict{Tuple{Node, Subgraph}, ProbabilityDistribution}(),
                                            Subgraph[])
end

"""
Return currently active RecognitionFactorization.
Create one if there is none.
"""
function currentRecognitionFactorization()
    try
        return current_recognition_factorization
    catch
        return RecognitionFactorization()
    end
end

setCurrentRecognitionFactorization(rf::RecognitionFactorization) = global current_recognition_factorization = rf

RecognitionFactorization() = setCurrentRecognitionFactorization(RecognitionFactorization(currentGraph()))

show(io::IO, rf::RecognitionFactorization) = println(io, "RecognitionFactorization with $(length(subgraphs)) factors")

function generateSubgraphId(rf::RecognitionFactorization=currentRecognitionFactorization())
    # Automatically generates a unique subgraph id
    rf.subgraph_id_counter += 1
    return Symbol("subgraph$(rf.subgraph_id_counter)")
end

# Subgraph getters
subgraph(id::Symbol, rf::RecognitionFactorization=currentRecognitionFactorization()) = rf.subgraphs[id]
subgraphs(ids::Vector{Symbol}, rf::RecognitionFactorization=currentRecognitionFactorization()) = Subgraph[subgraph(id, rf) for id in ids]
sg(id::Symbol, rf::RecognitionFactorization=currentRecognitionFactorization()) = subgraph(id, rf)
sg(ids::Vector{Symbol}, rf::RecognitionFactorization=currentRecognitionFactorization()) = subgraphs(ids, rf)

"""
Add a subgraph to the recognition factorization
A factor is contructed for each set/tuple (structured) or edge (naive) in the argument list
"""
function factor(edge_set::Set{Edge}, rf::RecognitionFactorization=currentRecognitionFactorization(); id::Symbol=generateSubgraphId(rf))
    # Construct a new Subgraph
    internal_edges = extend(edge_set)
    external_edges = setdiff(edges(nodes(internal_edges)), internal_edges) # external_edges are the difference between all edges connected to nodes, and the internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), nodes(internal_edges)) # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges

    subgraph = Subgraph(id, internal_edges, external_edges, sort(collect(nodes_connected_to_external_edges)))
    rf.subgraphs[id] = subgraph # Add subgraph to factorization
    
    # Initialize lookup tables
    for node in nodes_connected_to_external_edges
        rf.node_subgraph_to_internal_edges[(node, subgraph)] = intersect(edges(node), internal_edges)
    end
    for internal_edge in internal_edges
        rf.edge_to_subgraph[internal_edge] = subgraph
    end

    # Add the subgraph to the default update order
    push!(rf.default_update_order, subgraph)

    return subgraph
end

factor(edge::Edge, rf::RecognitionFactorization=currentRecognitionFactorization(); id::Symbol=generateSubgraphId(rf)) = factor(Set([edge]), rf; id=id)

# Tuple specifies a structured factor
factor(edge_tuple::Tuple, rf::RecognitionFactorization=currentRecognitionFactorization(); id::Symbol=generateSubgraphId(rf)) = factor(Set(edge_tuple), rf; id=id)

# Vector specifies a separate factor for each entry (Edge, tuple of Set{Edge})
function factor(edge_factors::Vector, rf::RecognitionFactorization=currentRecognitionFactorization(); ids::Vector{Symbol}=Symbol[])
    for (k, edge_factor) in enumerate(edge_factors)
        isempty(ids) ? id = generateSubgraphId(rf) : id = ids[k]
        factor(edge_factor, rf; id=id)
    end
end

"""
Find the smallest legal subgraph (connected through deterministic nodes) that includes the argument edges
"""
function extend(edge_set::Set{Edge})
    cluster = Set{Edge}() # Set to fill with edges in equality cluster
    edges = copy(edge_set)
    while length(edges) > 0 # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(cluster, current_edge) # Add to edge cluster
        for node in [current_edge.head.node, current_edge.tail.node] # Check both head and tail node for deterministic type
            if ForneyLab.isDeterministic(node)
                for interface in node.interfaces
                    if !is(interface.edge, current_edge) && !(interface.edge in cluster) # Is next level edge not seen yet?
                        push!(edges, interface.edge) # Add to buffer to visit sometime in the future
                    end
                end
            end
        end
    end

    return cluster
end

extend(edge::Edge) = extend(Set([edge]))

"""
Apply mean field factorization
"""
function factorizeMeanField(rf::RecognitionFactorization=currentRecognitionFactorization())
    unfactored_edges = copy(edges(rf.graph)) # Iterate through ordered vector of edges
    while length(unfactored_edges) > 0
        current_edge = sort(collect(unfactored_edges))[1] # Pick the next edge from an ordered set of edges; this ensures that the factorization is deterministic
        subgraph = factor(current_edge, rf)
        cluster = subgraph.internal_edges # Get edges that were just added to factor
        setdiff!(unfactored_edges, cluster) # Remove factored edges
    end

    return rf
end

"""
Define the initial recognition distribution for an Edge, Set{Edge}, tuple of edges or Node-Subgraph combination
"""
function initialize(node::Node, sg::Subgraph, initial_rd::ProbabilityDistribution, rf::RecognitionFactorization=currentRecognitionFactorization())
    # Set initial distribution
    rf.initial_recognition_distributions = Partitioned(push!(rf.initial_recognition_distributions.factors, initial_rd))
    
    # Set iteration distribution
    rf.recognition_distributions = Partitioned(push!(rf.recognition_distributions.factors, deepcopy(initial_rd)))
    rf.node_subgraph_to_recognition_distribution[(node, sg)] = deepcopy(initial_rd)

    return rf
end

function initialize(edges::Set{Edge}, initial_rd::ProbabilityDistribution, rf::RecognitionFactorization=currentRecognitionFactorization())
    # Find the subgraph that the edge set belongs to
    subgraphs = unique([rf.edge_to_subgraph[edge] for edge in edges])
    (length(subgraphs) == 0) && error("Cannot specify recognition distribution over unfactorized edge(s). Use factor to specify a recognition factorization.")
    (length(subgraphs) == 1) || error("Cannot specify recognition distribution when edges are distributed over multiple subgraphs")
    subgraph = subgraphs[1]

    # Find the external node in the subgraph that is connected to the edge set
    node_set = intersect(nodes(edges), Set(subgraph.nodes_connected_to_external_edges))
    (length(node_set) == 0) && error("Cannot specify recognition distribution over specified edge(s)")
    (length(node_set) > 2) && error("Cannot specify recognition distribution over specified edge(s)")
    for node in node_set # Note: an edge could be connected to two external nodes, hence the loop 
        initialize(node, subgraph, initial_rd, rf)
    end

    return rf
end

initialize(edge::Edge, initial_rd::ProbabilityDistribution, rf::RecognitionFactorization=currentRecognitionFactorization()) = initialize(Set([edge]), initial_rd, rf)

# Specify a joint recognition distribution
initialize(edges::Tuple, initial_rd::ProbabilityDistribution, rf::RecognitionFactorization=currentRecognitionFactorization()) = initialize(Set(edges), initial_rd, rf)

# A vector argument initializes a recognition distribution for each entry
function initialize(edge_clusters::Vector, initial_rd::ProbabilityDistribution, rf::RecognitionFactorization=currentRecognitionFactorization())
    for edge_cluster in edge_clusters
        initialize(edge_cluster, initial_rd, rf)
    end
end

"""
Verify whether the recognition factorization is proper and all distributions are set.
Throw an error if the factorization is improper.
"""
function verifyProper(rf::RecognitionFactorization)
    # All edges in the graph should be internal to a subgraph;
    # And no edge should be internal to multiple subgraphs
    checked_edges = Edge[]
    for edge in edges(rf.graph)
        haskey(rf.edge_to_subgraph, edge) || error("Edge $(edge) is not yet associated with a subgraph")
        (edge in checked_edges) && error("Edge $(edge) is associated with multiple subgraphs")
        push!(checked_edges, edge)
    end

    # All node-subgraph combinations should be coupled to a recognition distribution;
    # And all node-subgraph combinations should be coupled to a set of internal edges
    for subgraph in values(rf.subgraphs)
        for external_node in subgraph.nodes_connected_to_external_edges
            haskey(rf.node_subgraph_to_recognition_distribution, (external_node, subgraph)) || error("A node-subgraph combination for $(external_node) is not yet associated with a recognition distribution")
            haskey(rf.node_subgraph_to_internal_edges, (external_node, subgraph)) || error("A node-subgraph combination for $(external_node) is not yet associated with internal edges")
        end
    end

    return true
end

"""
Reset the recognition distributions to the initial distributions
"""
function resetRecognitionDistributions!(rf::RecognitionFactorization)
    n_factors = length(rf.recognition_distributions.factors)
    for factor_ind in 1:n_factors
        initial_rd = rf.initial_recognition_distributions.factors[factor_ind]
        rd = rf.recognition_distributions.factors[factor_ind]
        injectParameters!(rd, initial_rd) # Inject parametesr of initial distribution into recognition distribution
    end

    return rf
end