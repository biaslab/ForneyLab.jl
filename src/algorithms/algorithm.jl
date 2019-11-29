export Algorithm, currentAlgorithm

"""
An `Algorithm` holds a collection of (non-overlapping) recognition factors that
specify the recognition factorization over a factor graph that is used for variational inference.
"""
mutable struct Algorithm
    graph::FactorGraph
    recognition_factors::Dict{Symbol, RecognitionFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_recognition_factor::Dict{Edge, RecognitionFactor}
    node_edge_to_cluster::Dict{Tuple{FactorNode, Edge}, Cluster}

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

Algorithm() = setCurrentAlgorithm(
    Algorithm(
        currentGraph(),
        Dict{Symbol, RecognitionFactor}(),
        Dict{Edge, RecognitionFactor}(),
        Dict{Tuple{FactorNode, Edge}, Symbol}(),
        Dict{Symbol, Any}[],
        Dict{Symbol, Any}[]))

"""
Construct a `Algorithm` consisting of one
`RecognitionFactor` for each argument
"""
function Algorithm(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[])
    rf = Algorithm()
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

function targetToMarginalEntry(algo::Algorithm)
    mapping = Dict{Union{Cluster, Variable}, MarginalEntry}()
    for rf in algo.recognition_factors
        rf_mapping = targetToMarginalEntry(rf.marginal_table)
        merge!(mapping, rf_mapping)
    end

    return mapping    
end