export RecognitionFactor

"""
A `RecognitionFactor` specifies the subset of variables that comprise
a joint factor in the recognition factorization. A RecognitionFactor can
be defined by a (selection of) variable(s), or on the complete graph.
"""
mutable struct RecognitionFactor
    id::Symbol
    internal_edges::Set{Edge}

    # Fields set by algorithm construction
    target_variables::Set{Variable} # Target variables for which marginals are required
    target_clusters::Set{Cluster} # Target clusters for which marginals are required

    # Fields set by algorithm assembler
    algorithm_id::Symbol # Specify the algorithm id for this recognition_factor
    schedule::Schedule # Specify the internal message passing schedule for this recognition factor
    marginal_table::MarginalTable # Specify the marginal updates for internal variables
    optimize::Bool # Indicate the need for an optimization block
    initialize::Bool # Indicate the need for a message initialization block

    function RecognitionFactor(algo=currentAlgorithm(); id=generateId(RecognitionFactor))
        internal_edges = nonClampedEdges(algo.graph) # Include all non-clamped edges in a single recognition factor
        self = new(id, internal_edges)
        algo.recognition_factors[id] = self # Register self with the algorithm

        return self
    end

    function RecognitionFactor(variables::Set{Variable}; algo=currentAlgorithm(), id=generateId(RecognitionFactor))
        internal_edges = extend(edges(variables)) # Include all deterministically liked variables in a single recognition factor
        self = new(id, internal_edges)
        algo.recognition_factors[id] = self # Register self with the algorithm

        return self
    end
end

RecognitionFactor(variable::Variable; algo=currentAlgorithm(), id=generateId(RecognitionFactor)) = RecognitionFactor(Set([variable]), algo=algo, id=id)
RecognitionFactor(variables::Vector{Variable}; algo=currentAlgorithm(), id=generateId(RecognitionFactor)) = RecognitionFactor(Set(variables), algo=algo, id=id)

"""
Find edges that are internal to the recognition factor and connected to node.
This function is used for constructing clusters. Therefore, the edges are returned 
in the same order as the node's interfaces. Optionally ignores clamped edges.
"""
function localInternalEdges(node::FactorNode, rf::RecognitionFactor)
    local_internal_edges = Edge[]
    for interface in node.interfaces
        if interface.edge in rf.internal_edges # Edge is internal to rf
            push!(local_internal_edges, interface.edge) # Otherwise include the edge
        end
    end

    return local_internal_edges
end

function draw(rf::RecognitionFactor; schedule=ScheduleEntry[], args...)
    subgraph_nodes = nodes(rf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), rf.internal_edges)
    ForneyLab.graphviz(ForneyLab.genDot(subgraph_nodes, rf.internal_edges, schedule=schedule, external_edges=external_edges); args...)
end

"""
Find the smallest legal subgraph that includes the argument edges. Default setting terminates the search at soft factors
and does not constrain the search to a limiting set (as specified by an empty `limit_set` argument).
"""
function extend(edge_set::Set{Edge}; terminate_at_soft_factors=true, limit_set=Set{Edge}(), include_clamped=false)
    extension = Set{Edge}() # Initialize extension
    edges = copy(edge_set) # Initialize stack
    while !isempty(edges) # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(extension, current_edge) # Add to edge extension
        
        connected_nodes = [] # Find nodes connected to edge (as a vector)
        (current_edge.a == nothing) || push!(connected_nodes, current_edge.a.node)
        (current_edge.b == nothing) || push!(connected_nodes, current_edge.b.node)

        for node in connected_nodes # Check both head and tail node (if present)
            if (terminate_at_soft_factors==false) || isa(node, DeltaFactor)
                for interface in node.interfaces
                    partner = ultimatePartner(interface)
                    if !include_clamped && (partner != nothing) && isa(partner.node, Clamp) # If clamps should be excluded, and there is a partner, and the partner is clamped
                        continue # Skip clamped edge (if desired)
                    end
                    if (interface.edge !== current_edge) && !(interface.edge in extension) && ( isempty(limit_set) || (interface.edge in limit_set) ) # No backtracking, if edge is not already visited and edge is contained within limit set
                        push!(edges, interface.edge) # Add unseen edges to the stack (to visit sometime in the future)
                    end
                end
            end
        end
    end

    return extension
end

extend(edge::Edge; terminate_at_soft_factors=true, limit_set=Set{Edge}(), include_clamped=false) = extend(Set{Edge}([edge]), terminate_at_soft_factors=terminate_at_soft_factors, limit_set=limit_set, include_clamped=include_clamped)

"""
Find all edges that are not clamped
"""
function nonClampedEdges(graph::FactorGraph)
    clamped_edges = Set{Edge}()
    for (id, node) in graph.nodes # Find edges terminated by a Clamp
        if isa(node, Clamp)
            push!(clamped_edges, node.interfaces[1].edge)
        end
    end

    non_clamped_edges = setdiff(edges(graph), clamped_edges)

    return non_clamped_edges
end

"""
Find the `RecognitionFactor` that `edge` belongs to (if available)
"""
function recognitionFactor(edge::Edge)
    dict = current_algorithm.edge_to_recognition_factor
    if haskey(dict, edge)
        rf = dict[edge]
    else # No recognition factor is found, return the edge itself
        rf = edge
    end

    return rf::Union{RecognitionFactor, Edge}
end

"""
Return the ids of the recognition factors to which edges connected to `node` belong
"""
localRecognitionFactors(node::FactorNode) = Any[recognitionFactor(interface.edge) for interface in node.interfaces]

"""
Find the nodes in `recognition_factor` that are connected to external edges
"""
function nodesConnectedToExternalEdges(rf::RecognitionFactor)
    subgraph_nodes = nodes(rf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), rf.internal_edges)
    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)

    return nodes_connected_to_external_edges
end

"""
Return a dictionary from recognition factors to regions local to `node`
"""
function localRecognitionFactorToRegion(node::FactorNode)
    # For each edge connected to node, collect the respective recognition factors and regions
    local_recognition_factors = localRecognitionFactors(node)
    local_regions = localRegions(node)

    # Construct dictionary for local recognition factorization
    local_recognition_factor_to_region = Dict{Union{RecognitionFactor, Edge}, Region}()
    for (idx, factor) in enumerate(local_recognition_factors)
        local_recognition_factor_to_region[factor] = local_regions[idx]
    end

    return local_recognition_factor_to_region
end

function collectAverageEnergyInbounds(node::FactorNode)
    inbounds = Any[]
    local_recognition_factor_to_region = localRecognitionFactorToRegion(node)

    encountered_recognition_factors = Union{RecognitionFactor, Edge}[] # Keep track of encountered recognition factors
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        current_recognition_factor = recognitionFactor(node_interface.edge)

        if (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(copy(inbound_interface.node), ProbabilityDistribution)) # Copy Clamp before assembly to prevent overwriting dist_or_msg field
        elseif !(current_recognition_factor in encountered_recognition_factors)
            # Collect marginal entry from marginal dictionary (if marginal entry is not already accepted)
            target = local_recognition_factor_to_region[current_recognition_factor]
            push!(inbounds, current_algorithm.target_to_marginal_entry[target])
        end

        push!(encountered_recognition_factors, current_recognition_factor)
    end

    return inbounds
end