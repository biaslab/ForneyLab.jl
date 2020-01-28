import Base:iterate, values
export 
RecognitionFactorization, 
currentRecognitionFactorization, 
setCurrentRecognitionFactorization

mutable struct RecognitionFactorization
    recognition_factors::Dict{Symbol, RecognitionFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_recognition_factor::Dict{Edge, RecognitionFactor}
    node_edge_to_cluster::Dict{Tuple{FactorNode, Edge}, Cluster}
end

"""
Return currently active `RecognitionFactorization`.
Create one if there is none.
"""
function currentRecognitionFactorization()
    try
        return current_recognition_factorization
    catch
        return RecognitionFactorization()
    end
end

function setCurrentRecognitionFactorization(rfz::RecognitionFactorization)          global current_recognition_factorization = rfz
end

function RecognitionFactorization() 
    setCurrentRecognitionFactorization(
        RecognitionFactorization(
            Dict{Symbol, RecognitionFactor}(),
            Dict{Edge, RecognitionFactor}(),
            Dict{Tuple{FactorNode, Edge}, Symbol}(),
        )
    )
end

iterate(rfz::RecognitionFactorization) = iterate(rfz.recognition_factors)
iterate(rfz::RecognitionFactorization, state) = iterate(rfz.recognition_factors, state)

function values(rfz::RecognitionFactorization)
    return values(rfz.recognition_factors)
end

"""
Construct a `RecognitionFactorization` consisting of one
`RecognitionFactor` for each argument
"""
function RecognitionFactorization(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[])
    rfz = RecognitionFactorization()
    isempty(ids) || (length(ids) == length(args)) || error("Length of ids must match length of recognition factor arguments")
    for (i, arg) in enumerate(args)
        if isempty(ids)
            RecognitionFactor(arg, id=generateId(RecognitionFactor))
        else        
            RecognitionFactor(arg, id=ids[i])
        end
    end
    return rfz
end
