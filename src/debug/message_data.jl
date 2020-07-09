export MessageData

mutable struct MessageData{T <: Union{String, Symbol}}
    edgeID::T
    type::T
    marginal::ProbabilityDistribution

    MessageData{T}(edgeID::T, type::T, marginal::ProbabilityDistribution) where {T<:String} = new()
end
