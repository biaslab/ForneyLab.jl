export MarginalData, MarginalSnapshot

struct MarginalSnapshot
    variate::String
    family::String
    params::Dict
end

struct MarginalData
    id::Union{String, Symbol}
    edgeIDs::Vector{String}
    marginal::MarginalSnapshot
end

MarginalSnapshot(dist::ProbabilityDistribution{V, F}) where { F, V }     = MarginalSnapshot(string(V), string(F), dist.params)

function MarginalSnapshot(dist::ProbabilityDistribution{V, F}) where { F <: Function, V } 
    # error("Cannot dump function marginal [WIP]")
    return MarginalSnapshot(string(V), string(F), Dict())
end
