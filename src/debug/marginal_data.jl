export MarginalData

struct MarginalData
    id::Union{String, Symbol}
    edgesIDs::Vector{String}
    marginal::ProbabilityDistribution
end
