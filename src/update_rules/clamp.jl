mutable struct SPClamp{T<:VariateType} <: SumProductRule{Clamp{T}} end
outboundType{T<:VariateType}(::Type{SPClamp{T}}) = Message{PointMass, T}
isApplicable{T<:VariateType}(::Type{SPClamp{T}}, input_types::Vector{<:Type}) = true