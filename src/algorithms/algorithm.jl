export Algorithm, currentAlgorithm

"""
An `Algorithm` holds a collection of (non-overlapping) recognition factors that
specify the recognition factorization over a factor graph.
"""
mutable struct Algorithm
    id::Symbol
    
    graph::FactorGraph
    recognition_factorization::RecognitionFactorization

    # Bookkeeping for faster lookup during assembly
    interface_to_schedule_entry::Dict{Interface, ScheduleEntry}
    target_to_marginal_entry::Dict{Union{Variable, Cluster}, MarginalEntry}

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

setCurrentAlgorithm(algo::Algorithm) = global current_algorithm = algo

Algorithm(id=Symbol("")) = setCurrentAlgorithm(
    Algorithm(
        id,
        currentGraph(),
        currentRecognitionFactorization(),
        Dict{Tuple{FactorNode, Edge}, Symbol}(),
        Dict{Interface, ScheduleEntry}(),
        Dict{Union{Variable, Cluster}, MarginalEntry}(),
        Dict{Symbol, Any}[],
        Dict{Symbol, Any}[]))


function interfaceToScheduleEntry(algo::Algorithm)
    mapping = Dict{Interface, ScheduleEntry}()
    for (id, rf) in algo.recognition_factorization
        rf_mapping = interfaceToScheduleEntry(rf.schedule)
        merge!(mapping, rf_mapping)
    end

    return mapping
end

function targetToMarginalEntry(algo::Algorithm)
    mapping = Dict{Union{Cluster, Variable}, MarginalEntry}()
    for (id, rf) in algo.recognition_factorization
        rf_mapping = targetToMarginalEntry(rf.marginal_table)
        merge!(mapping, rf_mapping)
    end

    return mapping    
end