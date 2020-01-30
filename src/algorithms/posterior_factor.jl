export PosteriorFactor

"""
A `PosteriorFactor` specifies the subset of variables that comprise
a joint factor in the posterior factorization.
"""
mutable struct PosteriorFactor
    id::Symbol
    variables::Set{Variable}
    clusters::Set{Cluster}
    internal_edges::Set{Edge}

    # Fields set by algorithm assembler
    algorithm_id::Symbol # Specify the algorithm id for this posterior_factor
    schedule::Schedule # Specify the internal message passing schedule for this posterior factor
    marginal_table::MarginalTable # Specify the marginal updates for internal variables
    optimize::Bool # Indicate the need for an optimization block
    initialize::Bool # Indicate the need for a message initialization block

    function PosteriorFactor(pfz=currentPosteriorFactorization(); id=generateId(PosteriorFactor))
        # Constructor for empty container
        self = new(id)
        pfz.posterior_factors[id] = self # Register self with the algorithm

        return self
    end

    function PosteriorFactor(
        variables::Set{Variable};
        pfz=currentPosteriorFactorization(),
        id=generateId(PosteriorFactor))
        
        # Determine nodes connected to external edges
        internal_edges = ForneyLab.extend(edges(variables))
        subgraph_nodes = nodes(internal_edges)
        external_edges = setdiff(edges(subgraph_nodes), internal_edges)
        nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)
        
        # Determine variables required for variational updates
        internal_edges_connected_to_external_nodes = intersect(edges(nodes_connected_to_external_edges), internal_edges)
        posterior_variables = Set{Variable}([edge.variable for edge in internal_edges_connected_to_external_nodes])
        
        # Construct clusters required for (structured) variational updates
        posterior_clusters = Set{Cluster}()
        for node in nodes_connected_to_external_edges
            # Cluster edges must be ordered according to interfaces (therefore, intersect(edges(node), internal_edges) will not suffice)
            cluster_edges = Edge[]
            for interface in node.interfaces
                if interface.edge in internal_edges
                    push!(cluster_edges, interface.edge)
                end
            end

            if length(cluster_edges) > 1 # Constuct Cluster if multiple edges connected to node belong to the current subgraph
                push!(posterior_clusters, Cluster(node, cluster_edges))
            end
        end            

        # Create new posterior factor
        self = new(id, union(variables, posterior_variables), posterior_clusters, internal_edges)
        pfz.posterior_factors[id] = self # Register self with the algorithm

        # Register relevant edges with the algorithm for fast lookup during scheduling
        for edge in internal_edges_connected_to_external_nodes
            pfz.edge_to_posterior_factor[edge] = self
        end

        # Register clusters with the algorithm for fast lookup during scheduling
        for cluster in posterior_clusters
            for edge in cluster.edges
                pfz.node_edge_to_cluster[(cluster.node, edge)] = cluster
            end
        end 

        return self 
    end
end

PosteriorFactor(variable::Variable; id=generateId(PosteriorFactor)) = PosteriorFactor(Set([variable]), id=id)
PosteriorFactor(variables::Vector{Variable}; id=generateId(PosteriorFactor)) = PosteriorFactor(Set(variables), id=id)

function draw(pf::PosteriorFactor; schedule=ScheduleEntry[], args...)
    subgraph_nodes = nodes(pf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), pf.internal_edges)
    ForneyLab.graphviz(ForneyLab.genDot(subgraph_nodes, pf.internal_edges, schedule=schedule, external_edges=external_edges); args...)
end

"""
Return whether the subgraph contains a collider. If a collider is found, this will lead to conditional dependencies in the posterior distribution (posterior).
"""
function hasCollider(pf::PosteriorFactor)
    stack = copy(pf.internal_edges)
    while !isempty(stack)
        # Choose a maximal connected cluster in the subgraph
        seed_edge = first(stack)
        connected_cluster = extend(seed_edge, terminate_at_soft_factors=false, limit_set=pf.internal_edges) # Extend the seed edge to find a maximal connected cluster within the subgraph

        if hasCollider(connected_cluster)
            # If one of the connected clusters has a collider, the subgraph has a collider
            return true
        end

        stack = setdiff(stack, connected_cluster)
    end

    return false
end

"""
Return whether connected_cluster contains a collider. This function assumes the graph for connected_cluster is a connected tree.
"""
function hasCollider(connected_cluster::Set{Edge})
    nodes_connected_to_external_edges = nodesConnectedToExternalEdges(connected_cluster)

    # A prior node is a node for which only the outbound edge is internal to the cluster.
    # It represents a prior belief over a single variable in the present subgraph.
    n_prior_nodes = 0
    for node in nodes_connected_to_external_edges
        if node.interfaces[1].edge in connected_cluster # Check if the outbound edge is internal to the subgraph
            # If so, check if there are no other edges internal to the subgraph
            prior_flag = true
            for iface in node.interfaces
                (iface == node.interfaces[1]) && continue
                if iface.edge in connected_cluster # There is another edge besides the outbound in the cluster
                    prior_flag = false
                    break
                end
            end
            prior_flag && (n_prior_nodes += 1) # If there is no other edge internal to the subgraph, the node is a prior node
        end

        # If the subgraph contains more than one prior node, the variables on the outbound
        # edges of these nodes are conditionally dependent in the posterior distibution (posterior)
        (n_prior_nodes > 1) && return true
    end

    return false
end

"""
Find the smallest legal subgraph that includes the argument edges. Default setting terminates the search at soft factors
and does not constrain the search to a limiting set (as specified by an empty `limit_set` argument).
"""
function extend(edge_set::Set{Edge}; terminate_at_soft_factors=true, limit_set=Set{Edge}())
    cluster = Set{Edge}() # Set to fill with edges in cluster
    edges = copy(edge_set)
    while !isempty(edges) # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(cluster, current_edge) # Add to edge cluster
        
        connected_nodes = [] # Find nodes connected to edge (as a vector)
        (current_edge.a == nothing) || push!(connected_nodes, current_edge.a.node)
        (current_edge.b == nothing) || push!(connected_nodes, current_edge.b.node)

        for node in connected_nodes # Check both head and tail node (if present)
            if (terminate_at_soft_factors==false) || isa(node, DeltaFactor)
                for interface in node.interfaces
                    if (interface.edge !== current_edge) && !(interface.edge in cluster) && ( isempty(limit_set) || (interface.edge in limit_set) ) # Is next level edge not seen yet, and is it contained in the limiting set?
                        # Add unseen edges to the stack (to visit sometime in the future)
                        push!(edges, interface.edge)
                    end
                end
            end
        end
    end

    return cluster
end

extend(edge::Edge; terminate_at_soft_factors=true, limit_set=Set{Edge}()) = extend(Set{Edge}([edge]), terminate_at_soft_factors=terminate_at_soft_factors, limit_set=limit_set)

"""
Find the `PosteriorFactor` that `edge` belongs to (if available)
"""
function PosteriorFactor(edge::Edge)
    dict = current_posterior_factorization.edge_to_posterior_factor
    if haskey(dict, edge)
        pf = dict[edge]
    else # No posterior factor is found, return the edge itself
        pf = edge
    end

    return pf::Union{PosteriorFactor, Edge}
end

"""
Return the ids of the posterior factors to which edges connected to `node` belong
"""
localPosteriorFactors(node::FactorNode) = Any[PosteriorFactor(interface.edge) for interface in node.interfaces]

"""
Find the nodes in `posterior_factor` that are connected to external edges
"""
nodesConnectedToExternalEdges(posterior_factor::PosteriorFactor) = nodesConnectedToExternalEdges(posterior_factor.internal_edges)

"""
Find the nodes connected to `internal_edges` that are also connected to external edges
"""
function nodesConnectedToExternalEdges(internal_edges::Set{Edge})
    subgraph_nodes = nodes(internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), internal_edges)
    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)

    return nodes_connected_to_external_edges
end

"""
Return a dictionary from posterior factor-id to variable/cluster-ids local to `node`
"""
function localPosteriorFactorization(node::FactorNode)
    # For each edge connected to node, collect the posterior factor and cluster id
    local_posterior_factors = localPosteriorFactors(node)
    local_clusters = localClusters(node)

    # Construct dictionary for local posterior factorization
    local_posterior_factorization = Dict{Union{PosteriorFactor, Edge}, Union{Cluster, Variable}}()
    for (idx, factor) in enumerate(local_posterior_factors)
        local_posterior_factorization[factor] = local_clusters[idx]
    end

    return local_posterior_factorization
end