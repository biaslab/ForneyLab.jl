import Base:iterate, values
export 
PosteriorFactorization, 
currentPosteriorFactorization, 
setCurrentPosteriorFactorization

mutable struct PosteriorFactorization
    posterior_factors::Dict{Symbol, PosteriorFactor}

    # Bookkeeping for faster lookup during scheduling
    edge_to_posterior_factor::Dict{Edge, PosteriorFactor}
    node_edge_to_cluster::Dict{Tuple{FactorNode, Edge}, Cluster}
end

"""
Return currently active `PosteriorFactorization`.
Create one if there is none.
"""
function currentPosteriorFactorization()
    try
        return current_posterior_factorization
    catch
        return PosteriorFactorization()
    end
end

function setCurrentPosteriorFactorization(pfz::PosteriorFactorization)          global current_posterior_factorization = pfz
end

function PosteriorFactorization() 
    setCurrentPosteriorFactorization(
        PosteriorFactorization(
            Dict{Symbol, PosteriorFactor}(),
            Dict{Edge, PosteriorFactor}(),
            Dict{Tuple{FactorNode, Edge}, Symbol}(),
        )
    )
end

iterate(pfz::PosteriorFactorization) = iterate(pfz.posterior_factors)
iterate(pfz::PosteriorFactorization, state) = iterate(pfz.posterior_factors, state)

function values(pfz::PosteriorFactorization)
    return values(pfz.posterior_factors)
end

"""
Construct a `PosteriorFactorization` consisting of one
`PosteriorFactor` for each argument
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
