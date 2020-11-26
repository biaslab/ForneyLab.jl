export PosteriorFactor

"""
A `PosteriorFactor` specifies the subset of variables that comprise
a joint factor in the posterior factorization. A PosteriorFactor can
be defined by a (selection of) variable(s), or on the complete graph.
"""
mutable struct PosteriorFactor
    id::Symbol
    internal_edges::Set{Edge}

    # Fields set by algorithm construction
    target_variables::Set{Variable} # Target variables for which marginals are required
    target_clusters::Set{Cluster} # Target clusters for which marginals are required
    breaker_interfaces::Set{Interface} # Target interfaces for which messages are required
    ep_sites::Set{Interface} # Sites for expectation propagation update

    # Fields set by algorithm assembler
    algorithm_id::Symbol # Specify the algorithm id for this posterior_factor
    schedule::Schedule # Specify the internal message passing schedule for this posterior factor
    marginal_table::MarginalTable # Specify the marginal updates for internal variables
    optimize::Bool # Indicate the need for an optimization block
    initialize::Bool # Indicate the need for a message initialization block

    # Contruct a posterior factor that encompasses all variables in the graph
    function PosteriorFactor(fg::FactorGraph; pfz=currentPosteriorFactorization(), id=generateId(PosteriorFactor))
        internal_edges = edges(fg) # Include all edges in a single posterior factor
        self = new(id, internal_edges, Set{Variable}(), Set{Cluster}(), Set{Interface}(), Set{Interface}()) # Initialize posterior factor
        pfz.posterior_factors[id] = self # Register self with the algorithm

        return self
    end

    # Construct a posterior factor that includes all variables deterministically linked to the target variables
    function PosteriorFactor(seed_variables::Set{Variable}; pfz=currentPosteriorFactorization(), id=generateId(PosteriorFactor))
        internal_edges = extend(edges(seed_variables)) # Include all deterministically liked variables in a single posterior factor
        self = new(id, internal_edges, Set{Variable}(), Set{Cluster}(), Set{Interface}(), Set{Interface}()) # Initialize posterior factor
        pfz.posterior_factors[id] = self # Register self with the posterior factorization

        return self
    end
end

PosteriorFactor(seed_variable::Variable; pfz=currentPosteriorFactorization(), id=generateId(PosteriorFactor)) = PosteriorFactor(Set([seed_variable]), pfz=pfz, id=id)
PosteriorFactor(seed_variables::Vector{Variable}; pfz=currentPosteriorFactorization(), id=generateId(PosteriorFactor)) = PosteriorFactor(Set(seed_variables), pfz=pfz, id=id)

"""
`messagePassingSchedule()` generates a message passing schedule for the posterior factor
"""
function messagePassingSchedule(pf::PosteriorFactor)
    nodes_connected_to_external_edges = nodesConnectedToExternalEdges(pf)

    # Schedule messages towards posterior distributions and target sites, limited to the internal edges
    breaker_interfaces = sort(collect(pf.breaker_interfaces), rev=true)
    schedule = summaryPropagationSchedule(sort(collect(pf.target_variables), rev=true), 
                                          sort(collect(pf.target_clusters), rev=true),
                                          limit_set=pf.internal_edges,
                                          target_sites=breaker_interfaces)
    for entry in schedule
        if entry.interface in pf.ep_sites
            entry.message_update_rule = ExpectationPropagationRule{typeof(entry.interface.node)}
        elseif (entry.interface.node in nodes_connected_to_external_edges) && !isa(entry.interface.node, DeltaFactor)
            local_posterior_factors = localPosteriorFactors(entry.interface.node)
            if allunique(local_posterior_factors) # Local posterior factorization is naive
                entry.message_update_rule = NaiveVariationalRule{typeof(entry.interface.node)}
            else
                entry.message_update_rule = StructuredVariationalRule{typeof(entry.interface.node)}
            end
        else
            entry.message_update_rule = SumProductRule{typeof(entry.interface.node)}
        end
    end

    breaker_types = breakerTypes(collect(pf.breaker_interfaces))
    inferUpdateRules!(schedule, inferred_outbound_types=breaker_types)

    return schedule
end

function draw(pf::PosteriorFactor; schedule=ScheduleEntry[], args...)
    subgraph_nodes = nodes(pf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), pf.internal_edges)
    ForneyLab.graphviz(ForneyLab.genDot(subgraph_nodes, pf.internal_edges, schedule=schedule, external_edges=external_edges); args...)
end

"""
Find the smallest legal subgraph that includes the argument edges. Default setting terminates the search at soft factors
and does not constrain the search to a limiting set (as specified by an empty `limit_set` argument).
"""
function extend(edge_set::Set{Edge}; terminate_at_soft_factors=true, limit_set=Set{Edge}())
    extension = Set{Edge}() # Initialize extension
    stack = copy(edge_set) # Initialize stack
    while !isempty(stack) # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(stack) # Pick one
        push!(extension, current_edge) # Add to edge extension
        
        connected_nodes = [] # Find nodes connected to edge (as a vector)
        (current_edge.a == nothing) || push!(connected_nodes, current_edge.a.node)
        (current_edge.b == nothing) || push!(connected_nodes, current_edge.b.node)

        for node in connected_nodes # Check both head and tail node (if present)
            if !terminate_at_soft_factors || isa(node, DeltaFactor)
                for interface in node.interfaces
                    if (interface.edge !== current_edge) && !(interface.edge in extension) && ( isempty(limit_set) || (interface.edge in limit_set) ) # No backtracking, if edge is not already visited and edge is contained within limit set
                        push!(stack, interface.edge) # Add unseen edges to the stack (to visit sometime in the future)
                    end
                end
            end
        end
    end

    return extension
end

extend(edge::Edge; terminate_at_soft_factors=true, limit_set=Set{Edge}()) = extend(Set{Edge}([edge]), terminate_at_soft_factors=terminate_at_soft_factors, limit_set=limit_set)

"""
Find the `PosteriorFactor` that `edge` belongs to (if available)
"""
function posteriorFactor(edge::Edge)
    dict = current_posterior_factorization.edge_to_posterior_factor
    if haskey(dict, edge)
        pf = dict[edge]
    else # No posterior factor is found or edge is clamped; return the edge itself
        pf = edge
    end

    return pf::Union{PosteriorFactor, Edge}
end

"""
Check if edge or interface is terminated by a Clamp node
"""
isClamped(interface::Interface) = isdefined(interface, :node) && isa(interface.node, Clamp)
isClamped(edge::Edge) = isClamped(edge.a) || isClamped(edge.b)
isClamped(::Nothing) = false

"""
Return the ids of the posterior factors to which edges connected to `node` belong
"""
localPosteriorFactors(node::FactorNode) = Any[posteriorFactor(interface.edge) for interface in node.interfaces]

"""
Find the nodes in `posterior_factor` that are connected to external edges
"""
function nodesConnectedToExternalEdges(pf::PosteriorFactor)
    subgraph_nodes = nodes(pf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), pf.internal_edges)

    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)

    return nodes_connected_to_external_edges
end

"""
Return a dictionary from posterior factors to regions local to `node`
"""
function localEdgeToRegion(node::FactorNode)
    # For each edge connected to node, collect the respective posterior factors and regions
    local_edge_to_region = Dict{Edge, Region}()
    for interface in node.interfaces
        local_edge_to_region[interface.edge] = region(node, interface.edge)
    end

    return local_edge_to_region
end

function collectAverageEnergyInbounds(node::FactorNode)
    inbounds = Any[]
    local_edge_to_region = localEdgeToRegion(node)

    encountered_regions = Region[] # Keep track of encountered regions
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        current_region = region(node_interface.node, node_interface.edge)

        if isClamped(inbound_interface)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(copy(inbound_interface.node), ProbabilityDistribution)) # Copy Clamp before assembly to prevent overwriting dist_or_msg field
        elseif !(current_region in encountered_regions)
            # Collect marginal entry from marginal dictionary (if marginal entry is not already accepted)
            target = local_edge_to_region[node_interface.edge]
            push!(inbounds, current_inference_algorithm.target_to_marginal_entry[target])
        end

        push!(encountered_regions, current_region)
    end

    return inbounds
end