export RecognitionFactor, RecognitionFactorization, currentRecognitionFactorization

"""
A Cluster specifies a collection of `edges` adjacent to `node` that belong to the same
RecognitionFactor. A joint marginal can be computed over a cluster.
"""
mutable struct Cluster <: AbstractCluster
    id::Symbol
    node::FactorNode
    edges::Vector{Edge}

    function Cluster(node::FactorNode, edges::Vector{Edge})
        id = Symbol(join([edge.variable.id for edge in edges], "_"))
        self = new(id, node, edges)
        return self
    end
end

Base.isless(c1::Cluster, c2::Cluster) = isless("$(c1.id)", "$(c2.id)")

"""
A RecognitionFactor specifies the subset of variables that comprise
a joint factor in the recognition factorization.
"""
mutable struct RecognitionFactor
    id::Symbol
    variables::Set{Variable}
    clusters::Set{Cluster}
    internal_edges::Set{Edge}

    function RecognitionFactor(variables::Set{Variable}; rfz=currentRecognitionFactorization(), id=generateId(RecognitionFactor))
        # Determine nodes connected to external edges
        internal_edges = ForneyLab.extend(edges(variables))
        subgraph_nodes = nodes(internal_edges)
        external_edges = setdiff(edges(subgraph_nodes), internal_edges)
        nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)
        
        # Determine variables required for variational updates
        internal_edges_connected_to_external_nodes = intersect(edges(nodes_connected_to_external_edges), internal_edges)
        recognition_variables = Set{Variable}([edge.variable for edge in internal_edges_connected_to_external_nodes])
        
        # Construct clusters required for (structured) variational updates
        recognition_clusters = Set{Cluster}()
        for node in nodes_connected_to_external_edges
            # Cluster edges must be ordered according to interfaces (therefore, intersect(edges(node), internal_edges) will not suffice)
            cluster_edges = Edge[]
            for interface in node.interfaces
                if interface.edge in internal_edges
                    push!(cluster_edges, interface.edge)
                end
            end

            if length(cluster_edges) > 1 # Constuct Cluster if multiple edges connected to node belong to the current subgraph
                push!(recognition_clusters, Cluster(node, cluster_edges))
            end
        end            

        # Create new recognition factor
        self = new(id, union(variables, recognition_variables), recognition_clusters, internal_edges)
        rfz.recognition_factors[id] = self # Register self with recognition factorization

        # Register relevant edges with the recognition factorization for fast lookup during scheduling
        for edge in internal_edges_connected_to_external_nodes
            rfz.edge_to_recognition_factor[edge] = self
        end

        # Register clusters with the recognition factorization for fast lookup during scheduling
        for cluster in recognition_clusters
            for edge in cluster.edges
                rfz.node_edge_to_cluster[(cluster.node, edge)] = cluster
            end
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
                    if (interface.edge !== current_edge) && !(interface.edge in cluster) # Is next level edge not seen yet?
                        # If the node is a DeltaFactor, add unseen edges to the stack (to visit sometime in the future)
                        push!(edges, interface.edge)
                    end
                end
            end
        end
    end

    return cluster
end

"""
A RecognitionFactorization holds a collection of (non-overlapping) recognition factors that
specify the recognition factorization over a factor graph that is used for variational inference.
"""
mutable struct RecognitionFactorization
    graph::FactorGraph
    recognition_factors::Dict{Symbol, RecognitionFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_recognition_factor::Dict{Edge, RecognitionFactor}
    node_edge_to_cluster::Dict{Tuple{FactorNode, Edge}, Cluster}
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

RecognitionFactorization() = setCurrentRecognitionFactorization(
    RecognitionFactorization(
        currentGraph(),
        Dict{Symbol, RecognitionFactor}(),
        Dict{Edge, RecognitionFactor}(),
        Dict{Tuple{FactorNode, Edge}, Symbol}()))

"""
Construct a RecognitionFactorization consisting of one
RecognitionFactor for each argument
"""
function RecognitionFactorization(args...; ids=Symbol[])
    rf = RecognitionFactorization()
    isempty(ids) || (length(ids) == length(args)) || error("Length of ids must match length of recognition factor arguments")
    for (i, arg) in enumerate(args)
        if isempty(ids)
            RecognitionFactor(arg, id=generateId(RecognitionFactor))
        else        
            RecognitionFactor(arg, id=ids[i])
        end
    end
    return rf
end

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
Return the id of the cluster/variable that the node-edge combination belongs to 
"""
function clusterId(node::FactorNode, edge::Edge)
    dict = current_recognition_factorization.node_edge_to_cluster
    if haskey(dict, (node, edge))
        id = dict[(node, edge)].id
    else # No cluster is found, return the variable id
        id = edge.variable.id
    end

    return id
end

"""
Return the ids of the clusters/variables to which edges connected to `node` belong
"""
localClusterIds(node::FactorNode) = [clusterId(node, interface.edge) for interface in node.interfaces]

"""
Return a dictionary from recognition factor-id to variable/cluster-ids local to node
"""
function localRecognitionFactorization(node::FactorNode)
    # For each edge connected to node, collect the recognition factor and cluster id
    local_recognition_factor_ids = localRecognitionFactorIds(node)
    local_cluster_ids = localClusterIds(node)

    # Construct dictionary for local recognition factorization
    local_recognition_factorization = Dict{Symbol, Symbol}()
    for (idx, factor_id) in enumerate(local_recognition_factor_ids)
        local_recognition_factorization[factor_id] = local_cluster_ids[idx]
    end

    return local_recognition_factorization
end

"""
Find the nodes in `recognition_factor` that are connected to external edges
"""
function nodesConnectedToExternalEdges(recognition_factor::RecognitionFactor)
    internal_edges = recognition_factor.internal_edges
    subgraph_nodes = nodes(internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), internal_edges)
    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)

    return nodes_connected_to_external_edges
end