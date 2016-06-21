export RecognitionFactorization, currentRecognitionFactorization
export show, addFactor, factorizeMeanField, setRecognitionDistribution, setDefaultRecognitionDistributions

type RecognitionFactorization
    graph::FactorGraph

    # Recognition factorization
    subgraphs::Vector{Subgraph}
    edge_to_subgraph::Dict{Edge, Subgraph}
    node_subgraph_to_internal_edges::Dict{Tuple{Node, Subgraph}, Set{Edge}}

    # Recognition distributions
    recognition_distributions::Partitioned
    node_subgraph_to_recognition_distribution::Dict{Tuple{Node, Subgraph}, ProbabilityDistribution}

    RecognitionFactorization(graph) = new(  graph,
                                            Subgraph[],
                                            Dict{Edge, Subgraph}(),
                                            Dict{Tuple{Node, Subgraph}, Set{Edge}}(),
                                            Partitioned{ProbabilityDistribution,0}(ProbabilityDistribution[]),
                                            Dict{Tuple{Node, Subgraph}, ProbabilityDistribution}())
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

show(io::IO, rf::RecognitionFactorization) = println(io, "RecognitionFactorization with $(length(rf.subgraphs)) factors")

"""
Add a subgraph to the recognition factorization
"""
function addFactor(edge_set::Set{Edge}, rf::RecognitionFactorization=currentRecognitionFactorization())
    # Construct a new Subgraph
    internal_edges = extend(edge_set)
    external_edges = setdiff(edges(nodes(internal_edges)), internal_edges) # external_edges are the difference between all edges connected to nodes, and the internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), nodes(internal_edges)) # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges

    subgraph = Subgraph(internal_edges, external_edges, sort(collect(nodes_connected_to_external_edges)))
    push!(rf.subgraphs, subgraph) # Add subgraph to factorization
    
    # Initialize lookup tables
    for node in nodes_connected_to_external_edges
        rf.node_subgraph_to_internal_edges[(node, subgraph)] = intersect(edges(node), internal_edges)
    end
    for internal_edge in internal_edges
        rf.edge_to_subgraph[internal_edge] = subgraph
    end

    return rf
end

addFactor(edge::Edge, rf::RecognitionFactorization=currentRecognitionFactorization()) = addFactor(Set([edge]), rf)

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
        addFactor(current_edge, rf)
        cluster = rf.subgraphs[end].internal_edges # Get edges that were just added to factor
        setdiff!(unfactored_edges, cluster) # Remove factored edges
    end

    return rf
end

"""
Define the recognition distribution type for an Edge, Set{Edge} or Node-Subgraph combination
"""
function setRecognitionDistribution{T<:ProbabilityDistribution}(node::Node, sg::Subgraph, dist_type::Type{T}, rf::RecognitionFactorization=currentRecognitionFactorization())
    rd = vague(dist_type)
    rf.recognition_distributions = Partitioned(push!(rf.recognition_distributions.factors, rd)) # Adds a partition
    rf.node_subgraph_to_recognition_distribution[(node, sg)] = rd

    return rf
end

function setRecognitionDistribution{T<:ProbabilityDistribution}(edges::Set{Edge}, dist_type::Type{T}, rf::RecognitionFactorization=currentRecognitionFactorization())
    # Find the subgraph that the edge set belongs to
    subgraphs = unique([rf.edge_to_subgraph[edge] for edge in edges])
    (length(subgraphs) == 0) && error("Cannot specify recognition distribution over unfactorized edge(s). Use addFactor to specify a recognition factorization.")
    (length(subgraphs) == 1) || error("Cannot specify recognition distribution when edges are distributed over multiple subgraphs")
    subgraph = subgraphs[1]

    # Find the external node in the subgraph that is connected to the edge set
    node_set = intersect(nodes(edges), Set(subgraph.nodes_connected_to_external_edges))
    (length(node_set) == 0) && error("Cannot specify recognition distribution over specified edge(s)")
    (length(node_set) > 2) && error("Cannot specify recognition distribution over specified edge(s)")
    for node in node_set # Note: an edge could be connected to two external nodes, hence the loop 
        setRecognitionDistribution(node, subgraph, dist_type, rf)
    end

    return rf
end

setRecognitionDistribution{T<:ProbabilityDistribution}(edge::Edge, dist_type::Type{T}, rf::RecognitionFactorization=currentRecognitionFactorization()) = setRecognitionDistribution(Set([edge]), dist_type, rf)

# """
# Add default recognition distributions to the factorization
# """
# function setDefaultRecognitionDistributions(rf::RecognitionFactorization=currentRecognitionFactorization())
#     # TODO

#     return rf
# end

"""
Verify whether the recognition factorization proper and all distributions are set.
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
    for subgraph in rf.subgraphs
        for external_node in subgraph.nodes_connected_to_external_edges
            haskey(rf.node_subgraph_to_recognition_distribution, (external_node, subgraph)) || error("A node-subgraph combination for $(external_node) is not yet associated with a recognition distribution")
            haskey(rf.node_subgraph_to_internal_edges, (external_node, subgraph)) || error("A node-subgraph combination for $(external_node) is not yet associated with internal edges")
        end
    end

    return true
end

"""
Reset the recognition distributions to vague
"""
function resetRecognitionDistributions!(rf::RecognitionFactorization)
    for rd in rf.recognition_distributions.factors
        vague!(rd)
    end

    return rf
end

"""
Calculate the recognition distribution for the node-subgraph combination
"""
function calculateRecognitionDistribution!(rf::RecognitionFactorization, node::Node, subgraph::Subgraph)
    recognition_distribution = rf.node_subgraph_to_recognition_distribution[(node, subgraph)]
    internal_edges = rf.node_subgraph_to_internal_edges[(node, subgraph)]

    if length(internal_edges) == 1
        # Update for univariate q; when there is only one internal edge
        internal_edge = first(internal_edges)
        return prod!(internal_edge.tail.message.payload, internal_edge.head.message.payload, recognition_distribution)
    else
        # Update for multivariate q
        required_inputs = Array(Any, 0)
        for interface in node.interfaces # Iterate over all interfaces connected to node
            neighbouring_subgraph = rf.edge_to_subgraph[interface.edge]
            if neighbouring_subgraph == subgraph # edge is internal
                push!(required_inputs, interface.partner.message)
            else # edge is external
                push!(required_inputs, rf.node_subgraph_to_recognition_distribution[(node, neighbouring_subgraph)])
            end
        end
        return recognitionRule!(node, recognition_distribution, required_inputs...)
    end
end


############################
# Update rules for multivariate recognition distributions
############################

"""
NormalGamma update rule for variational message passing
"""
function recognitionRule!(  node::GaussianNode,
                            recognition_dist::NormalGamma,
                            msg_mean::Message{Gaussian},
                            msg_prec::Message{Gamma},
                            dist_y::Gaussian)

    dist_mean = msg_mean.payload
    dist_prec = msg_prec.payload
    ForneyLab.ensureParameters!(dist_mean, (:m,))
    ForneyLab.ensureParameters!(dist_y, (:W,))

    recognition_dist.m = dist_mean.m
    recognition_dist.beta = huge
    recognition_dist.a = dist_prec.a + 0.5
    recognition_dist.b = (1.0/(2.0*dist_y.W)) + dist_prec.b

    return recognition_dist
end