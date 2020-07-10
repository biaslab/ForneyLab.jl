export MessageData

struct MessageData
    edgeID::Union{String, Symbol}
    type::String
    marginal::ProbabilityDistribution
end
