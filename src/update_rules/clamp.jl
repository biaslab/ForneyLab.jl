mutable struct SPClamp{T<:VariateType} <: SumProductRule{Clamp{T}} end
outboundType(::Type{SPClamp{T}}) where T<:VariateType = Message{PointMass, T}
isApplicable(::Type{SPClamp{T}}, input_types::Vector{<:Type}) where T<:VariateType = true