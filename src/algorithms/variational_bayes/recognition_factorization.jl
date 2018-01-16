import Base.factor
export RecognitionFactor, RecognitionFactorization, currentRecognitionFactorization

"""
A RecognitionFactor specifies the subset of variables that comprise
a joint factor in the recognition factorization.
"""
type RecognitionFactor
    id::Symbol
    variables::Set{Variable}
    clusters::Set{Cluster}
    internal_edges::Set{Edge}

    function RecognitionFactor(variables::Set{Variable}; rfz=currentRecognitionFactorization(), id=generateId(RecognitionFactor))
        # Collect variables on internal edges connected to external nodes
        internal_edges = ForneyLab.extend(edges(variables))
        subgraph_nodes = nodes(internal_edges)
        external_edges = setdiff(edges(subgraph_nodes), internal_edges)
        nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)
        internal_edges_connected_to_external_nodes = intersect(edges(nodes_connected_to_external_edges), internal_edges)
        recognition_variables = Set{Variable}([edge.variable for edge in internal_edges_connected_to_external_nodes])
        clusters = Set{Cluster}(...) # TODO: determine clusters

        self = new(id, clusters, union(variables, recognition_variables), internal_edges)
        rfz.recognition_factors[id] = self # Register new factor with recognition factorization

        # Register internal edges with the recognition factorization for fast lookup during scheduling
        for edge in internal_edges_connected_to_external_nodes
            rfz.edge_to_recognition_factor[edge] = self
        end

        return self 
    end
end

RecognitionFactor(variable::Variable; id=generateId(RecognitionFactor)) = RecognitionFactor(Set([variable]), id=id)
RecognitionFactor(variables::Vector{Variable}; id=generateId(RecognitionFactor)) = RecognitionFactor(Set(variables), id=id)

function draw(rf::RecognitionFactor; schedule=ScheduleEntry[], args...)
    subgraph_nodes = nodes(rf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), rf.internal_edges)
    ForneyLab.graphviz(ForneyLab.genDot(subgraph_nodes, rf.internal_edges, schedule=schedule, external_edges=external_edges); args...)
end

"""
Find the smallest legal subgraph (connected through deterministic nodes) that includes the argument edges
"""
function extend(edge_set::Set{Edge})
    cluster = Set{Edge}() # Set to fill with edges in cluster
    edges = copy(edge_set)
    while length(edges) > 0 # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(cluster, current_edge) # Add to edge cluster
        for node in [current_edge.a.node, current_edge.b.node] # Check both head and tail node for deterministic type
            if isa(node, DeltaFactor)
                for interface in node.interfaces
                    if !is(interface.edge, current_edge) && !(interface.edge in cluster) # Is next level edge not seen yet?
                        # If the node is a DeltaFactor, add unseen edges to the stack (to visit sometime in the future)
                        push!(edges, interface.edge)
                    end
                end
            end
        end
    end

    return cluster
end

type RecognitionFactorization
    graph::FactorGraph
    recognition_factors::Dict{Symbol, RecognitionFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_recognition_factor::Dict{Edge, RecognitionFactor}
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

RecognitionFactorization() = setCurrentRecognitionFactorization(RecognitionFactorization(currentGraph(), Dict{Symbol, RecognitionFactor}(), Dict{Edge, RecognitionFactor}()))

"""
Return the id of the RecognitionFactor that `edge` belongs to
"""
function recognitionFactorId(edge::Edge)
    dict = current_recognition_factorization.edge_to_recognition_factor
    if haskey(dict, edge)
        rf_id = dict[edge].id
    else # No recognition factor is found, return a unique id
        rf_id = Symbol(name(edge))
    end

    return rf_id
end

"""
Return the ids of the recognition factors to which edges connected to `node` belong
"""
localRecognitionFactorIds(node::FactorNode) = [recognitionFactorId(interface.edge) for interface in node.interfaces]

"""
Return the ids of the clusters/variables to which edges connected to `node` belong
"""
function localClusterIds(node::FactorNode)
    # Note, a single edge might belong to multiple clusters
    for interface in node.interfaces
        ... # TODO
    end
end

"""
Return a dictionary from recognition factor-id to variable/cluster-ids local to node
"""
function localRecognitionFactorization(node::FactorNode)
    local_recognition_factor_ids = localRecognitionFactorIds(node)
    local_cluster_ids = localClusterIds(node)
    (length(local_recognition_factor_ids) == length(local_cluster_ids)) || error("Lengths of local recognition factorization and local clusters must agree")

    local_recognition_factorization = Dict{Symbol, Vector}()
    for (idx, factor) in enumerate(local_recognition_factor_ids)
        local_recognition_factorization[factor] = [local_cluster_ids[idx]]
    end

    return local_recognition_factorization
end

function nodesConnectedToExternalEdges(recognition_factor::RecognitionFactor)
    internal_edges = recognition_factor.internal_edges
    subgraph_nodes = nodes(internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), internal_edges)
    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)

    return nodes_connected_to_external_edges
end