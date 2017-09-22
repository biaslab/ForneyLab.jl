type SPClamp <: SumProductRule{Clamp} end
outboundType(::Type{SPClamp}) = Message{PointMass}
isApplicable(::Type{SPClamp}, input_types::Vector{DataType}) = true