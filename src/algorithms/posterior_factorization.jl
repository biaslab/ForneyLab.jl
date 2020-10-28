import Base:iterate, values
export 
PosteriorFactorization, 
currentPosteriorFactorization, 
setCurrentPosteriorFactorization

global current_posterior_factorization = nothing

mutable struct PosteriorFactorization
    graph::FactorGraph
    posterior_factors::Dict{Symbol, PosteriorFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_posterior_factor::Dict{Edge, PosteriorFactor}
    node_edge_to_cluster::Dict{Tuple{FactorNode, Edge}, Cluster}
    deterministic_edges::Set{Edge}

    # Fields for assembly
    energy_counting_numbers::Dict{FactorNode, Int64}
    entropy_counting_numbers::Dict{Region, Int64}

    # Flag to ensure required quantities for FE evaluation are computed
    free_energy_flag::Bool
end

"""
Return currently active `PosteriorFactorization`. Create one if there is none.
"""
function currentPosteriorFactorization()
    current_posterior_factorization === nothing ? PosteriorFactorization() : current_posterior_factorization
end

function setCurrentPosteriorFactorization(pfz::PosteriorFactorization)
    global current_posterior_factorization = pfz
end

"""
Initialize an empty `PosteriorFactorization` for sequential construction
"""
function PosteriorFactorization()
    setCurrentPosteriorFactorization(
        PosteriorFactorization(
            currentGraph(),
            Dict{Symbol, PosteriorFactor}(),
            Dict{Edge, PosteriorFactor}(),
            Dict{Tuple{FactorNode, Edge}, Symbol}(),
            deterministicEdges(currentGraph()),
            Dict{FactorNode, Int64}(),
            Dict{Region, Int64}(),
            false
        )
    )
end

"""
Initialize a `PosteriorFactorization` consisting of a single `PosteriorFactor` for the entire graph
"""
function PosteriorFactorization(fg::FactorGraph)
    pfz = setCurrentPosteriorFactorization(
        PosteriorFactorization(
            fg,
            Dict{Symbol, PosteriorFactor}(),
            Dict{Edge, PosteriorFactor}(),
            Dict{Tuple{FactorNode, Edge}, Symbol}(),
            deterministicEdges(fg),
            Dict{FactorNode, Int64}(),
            Dict{Region, Int64}(),
            false
        )
    )
    PosteriorFactor(fg, pfz=pfz, id=Symbol(""))
    return pfz
end

"""
Construct a `PosteriorFactorization` consisting of one `PosteriorFactor` for each argument
"""
function PosteriorFactorization(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[])
    pfz = PosteriorFactorization()
    isempty(ids) || (length(ids) == length(args)) || error("Length of ids must match length of posterior factor arguments")
    for (i, arg) in enumerate(args)
        if isempty(ids)
            PosteriorFactor(arg, id=generateId(PosteriorFactor))
        else        
            PosteriorFactor(arg, id=ids[i])
        end
    end
    return pfz
end

iterate(pfz::PosteriorFactorization) = iterate(pfz.posterior_factors)
iterate(pfz::PosteriorFactorization, state) = iterate(pfz.posterior_factors, state)
values(pfz::PosteriorFactorization) = values(pfz.posterior_factors)

"""
Populate the target regions of the PosteriorFactor.
The targets depend on the variables of interest, the local posterior
factorization, and whether free energy will be evaluated. At the same
time, fields for fast lookup during scheduling are populated in the
posterior factorization.
"""
function setTargets!(pf::PosteriorFactor, pfz::PosteriorFactorization; target_variables=Set{Variable}(), free_energy=false, external_targets=false)
    # Initialize empty set of target cluster node and edges.
    # We cannot build a Set of Clusters directly, because duplicate Clusters are not removed.
    large_regions = Set{Tuple}()

    # Register quantities of interest with the current posterior factor
    for variable in target_variables
        edge = first(variable.edges)
        if edge in pf.internal_edges
            push!(pf.target_variables, variable)
        end
    end

    # Determine which target regions are required by external posterior factors
    if external_targets
        nodes_connected_to_external_edges = nodesConnectedToExternalEdges(pf)
        for node in nodes_connected_to_external_edges
            isa(node, DeltaFactor) && continue # Skip deterministic nodes

            target_edges = localStochasticInternalEdges(node, pf, pfz) # Find internal edges connected to node
            if length(target_edges) == 1 # Only one internal edge, the marginal for a single Variable is required
                push!(pf.target_variables, target_edges[1].variable)
            elseif length(target_edges) > 1 # Multiple internal edges, register the region for computing the joint marginal
                push!(large_regions, (node, target_edges))
            end
        end
    end

    # Determine which targets are required for evaluating the free energy
    if free_energy
        # Initialize counting numbers
        variable_counting_numbers = Dict{Variable, Number}()
        cluster_counting_numbers = Dict{Tuple, Number}()

        # Iterate over large regions
        nodes_connected_to_internal_edges = nodes(pf.internal_edges)
        for node in nodes_connected_to_internal_edges
            target_edges = localStochasticInternalEdges(node, pf, pfz) # Find stochastic internal edges connected to node
            if isa(node, Clamp)
                continue
            elseif !isa(node, DeltaFactor) # Node is stochastic
                if length(target_edges) == 1 # Single internal edge
                    increase!(variable_counting_numbers, target_edges[1].variable, Inf) # For average energy evaluation, make sure to include the edge variable
                elseif length(target_edges) > 1 # Multiple internal edges
                    region = (node, target_edges) # Node and edges for region that will later define the Cluster
                    increase!(cluster_counting_numbers, region, Inf) # Make sure to include the region for average energy evaluation
                end
            elseif isa(node, Equality)
                increase!(variable_counting_numbers, target_edges[1].variable, 1) # Increase the counting number for the equality-constrained variable
            else # Node is deterministic and not Equality, requires a joint belief for the inbound variables
                outbound_edge = node.interfaces[1].edge
                inbound_target_edges = filter(edge->edge!=outbound_edge, target_edges) # Collect local internal stochastic inbound edges
                if length(inbound_target_edges) == 1 # Just a single local stochastic inbound edge
                    increase!(variable_counting_numbers, inbound_target_edges[1].variable, 1) # Increase counting number for variable on inbound edge
                elseif length(inbound_target_edges) >= 2 # Multiple stochastic inbound edges require the formation of a large (joint) region
                    region = (node, inbound_target_edges) # Region for node and stochastic inbound edges
                    increase!(cluster_counting_numbers, region, 1) # Increase the counting number for that region
                end
            end
        end

        # Iterate over small regions
        for edge in pf.internal_edges
            increase!(variable_counting_numbers, edge.variable, -(degree(edge) - 1)) # Discount single variables
        end

        # All targets with a non-zero counting number are required for free energy evaluation
        for (variable, cnt) in variable_counting_numbers
            if cnt != 0
                push!(pf.target_variables, variable)
            end
        end
        for (region, cnt) in cluster_counting_numbers
            if cnt != 0
                push!(large_regions, region)
            end
        end
    end

    # Determine which interfaces require breakers
    for edge in pf.internal_edges
        requiresBreaker(edge.a) && push!(pf.breaker_interfaces, edge.a)
        requiresBreaker(edge.b) && push!(pf.breaker_interfaces, edge.b)

        # Register internal edges with the posterior factorization for fast lookup during scheduling
        pfz.edge_to_posterior_factor[edge] = pf
    end
    
    # Create clusters, and register clusters with the posterior factorization for fast lookup during scheduling
    for region in large_regions
        cluster = Cluster(region...) # Use the region definition to construct a Cluster
        push!(pf.target_clusters, cluster)
        for edge in cluster.edges
            pfz.node_edge_to_cluster[(cluster.node, edge)] = cluster
        end
    end 

    # Register whether posterior factorization is prepared for free energy evaluation
    pfz.free_energy_flag = free_energy

    return pf
end

"""
Return the local stochastic regions around `node`
"""
function localStochasticRegions(node::FactorNode, pfz::PosteriorFactorization)
    regions = Region[]
    for interface in node.interfaces
        partner = ultimatePartner(interface)
        if !(interface.edge in pfz.deterministic_edges) # If edge is stochastic
            push!(regions, region(node, interface.edge))
        end
    end

    return regions
end

"""
Find stochastic edges that are internal to the posterior factor and connected to node.
This function is used for constructing clusters. Therefore, the edges are returned 
in the same order as the node's interfaces.
"""
function localStochasticInternalEdges(node::FactorNode, pf::PosteriorFactor, pfz::PosteriorFactorization)
    local_internal_edges = Edge[]
    for interface in node.interfaces
        if (interface.edge in pf.internal_edges) && !(interface.edge in pfz.deterministic_edges) # Edge is internal to posterior factor and stochastic
            push!(local_internal_edges, interface.edge) # Otherwise include the edge
        end
    end

    return local_internal_edges
end