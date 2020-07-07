mutable struct MessageData{T <: String}
    edgeID::T
    type::T
    marginal::ProbabilityDistribution

    MessageData{T}(edgeID::T, type::T, marginal::ProbabilityDistribution) where {T<:String} = new()
end
