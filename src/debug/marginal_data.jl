export MarginalData

mutable struct MarginalData
    edgesIDs::Vector{String}
    marginal::ProbabilityDistribution

    function MarginalData(edges, marginal)
        self = new(edges, marginal)
        return self
    end
end
