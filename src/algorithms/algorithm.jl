export Algorithm, currentAlgorithm

"""
An `Algorithm` holds a collection of (non-overlapping) recognition factors that
specify the recognition factorization over a factor graph.
"""
mutable struct Algorithm
    id::Symbol
    graph::FactorGraph
    recognition_factors::Dict{Symbol, RecognitionFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_recognition_factor::Dict{Edge, RecognitionFactor}
    node_edge_to_cluster::Dict{Tuple{FactorNode, Edge}, Cluster}

    # Bookkeeping for faster lookup during assembly
    interface_to_schedule_entry::Dict{Interface, ScheduleEntry}
    target_to_marginal_entry::Dict{Region, MarginalEntry}
    energy_counting_numbers::Dict{FactorNode, Int64}
    entropy_counting_numbers::Dict{Region, Int64}

    # Flag to ensure required quantities for FE evaluation are computed
    free_energy_flag::Bool

    # Fields for free energy algorithm assembly
    average_energies::Vector{Dict{Symbol, Any}}
    entropies::Vector{Dict{Symbol, Any}}
end

"""
Return currently active `Algorithm`.
Create one if there is none.
"""
function currentAlgorithm()
    try
        return current_algorithm
    catch
        return Algorithm()
    end
end

setCurrentAlgorithm(rf::Algorithm) = global current_algorithm = rf

Algorithm(id=Symbol("")) = setCurrentAlgorithm(
    Algorithm(
        id,
        currentGraph(),
        Dict{Symbol, RecognitionFactor}(),
        Dict{Edge, RecognitionFactor}(),
        Dict{Tuple{FactorNode, Edge}, Symbol}(),
        Dict{Interface, ScheduleEntry}(),
        Dict{Region, MarginalEntry}(),
        Dict{FactorNode, Int64}(),
        Dict{Region, Int64}(),
        false,
        Dict{Symbol, Any}[],
        Dict{Symbol, Any}[]))

"""
Construct a `Algorithm` consisting of one
`RecognitionFactor` for each argument
"""
function Algorithm(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[], id=Symbol(""))
    rf = Algorithm(id)
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
Pupulate the target regions of the RecognitionFactor.
The targets depend on the variables of interest, the local recognition
factorization, and whether free energy will be evaluated. At the same
time, fields for fast lookup during scheduling are populated in the algorithm.
"""
function setTargets!(rf::RecognitionFactor, algo::Algorithm, variables::Vector{Variable}=Variable[]; free_energy=false, external_targets=false)
    target_variables = Set{Variable}(variables) # Marginals of the quantities of interest are always required
    large_regions = Set{Tuple}() # Initialize empty set of target cluster node and edges. We cannot build a Set of Clusters directly, because duplicate Clusters are not removed.

    # Determine which target regions are required by external recognition factors
    if external_targets
        nodes_connected_to_external_edges = nodesConnectedToExternalEdges(rf)
        for node in nodes_connected_to_external_edges
            isa(node, DeltaFactor) && continue # Skip deterministic nodes

            target_edges = localInternalEdges(node, rf) # Find internal edges connected to node
            if length(target_edges) == 1 # Only one internal edge, the marginal for a single Variable is required
                push!(target_variables, target_edges[1].variable)
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
        nodes_connected_to_internal_edges = nodes(rf.internal_edges)
        for node in nodes_connected_to_internal_edges
            target_edges = localInternalEdges(node, rf) # Find internal edges connected to node
            if !isa(node, DeltaFactor) # Node is stochastic
                if length(target_edges) == 1 # Single internal edge
                    increase!(variable_counting_numbers, target_edges[1].variable, Inf) # For average energy evaluation, make sure to include the edge variable
                elseif length(target_edges) > 1 # Multiple internal edges
                    region = (node, target_edges) # Node and edges for region that will later define the Cluster,
                    increase!(cluster_counting_numbers, region, Inf) # Make sure to include the region for average energy evaluation
                end
            elseif isa(node, Equality)
                increase!(variable_counting_numbers, target_edges[1].variable, 1) # Increase the counting number for the equality-constrained variable
            else # Node is deterministic and not Equality
                if length(target_edges) == 2
                    increase!(variable_counting_numbers, target_edges[2].variable, 1) # Increase counting number for variable on inbound edge
                elseif length(target_edges) > 2
                    region = (node, target_edges[2:end]) # Region for node and inbound edges,
                    increase!(cluster_counting_numbers, region, 1) # Increase the counting number for that region
                end
            end
        end

        # Iterate over small regions
        for edge in rf.internal_edges
            increase!(variable_counting_numbers, edge.variable, -1) # Discount single variables
        end

        # All targets with a non-zero counting number are required for free energy evaluation
        for (variable, cnt) in variable_counting_numbers
            if cnt != 0
                push!(target_variables, variable)
            end
        end
        for (region, cnt) in cluster_counting_numbers
            if cnt != 0
                push!(large_regions, region)
            end
        end
    end

    # Register internal edges with the algorithm for fast lookup during scheduling
    for edge in rf.internal_edges
        algo.edge_to_recognition_factor[edge] = rf
    end

    # Create clusters, and register clusters with the algorithm for fast lookup during scheduling
    target_clusters = Set{Cluster}() # Initialize empty set of target clusters
    for region in large_regions # For each stored region definition
        cluster = Cluster(region...) # Use the region definition to construct a Cluster
        push!(target_clusters, cluster)
        for edge in cluster.edges
            algo.node_edge_to_cluster[(cluster.node, edge)] = cluster
        end
    end 

    # Register whether algorithm is prepared for free energy evaluation
    algo.free_energy_flag = free_energy

    # Register the targets with the recognition factor
    rf.target_variables = target_variables
    rf.target_clusters = target_clusters

    return rf
end

function interfaceToScheduleEntry(algo::Algorithm)
    mapping = Dict{Interface, ScheduleEntry}()
    for (id, rf) in algo.recognition_factors
        rf_mapping = interfaceToScheduleEntry(rf.schedule)
        merge!(mapping, rf_mapping)
    end

    return mapping
end

function targetToMarginalEntry(algo::Algorithm)
    mapping = Dict{Region, MarginalEntry}()
    for (id, rf) in algo.recognition_factors
        rf_mapping = targetToMarginalEntry(rf.marginal_table)
        merge!(mapping, rf_mapping)
    end

    return mapping    
end