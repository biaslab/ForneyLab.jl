type SPClamp <: SumProductRule{Clamp} end
outboundType(::Type{SPClamp}) = Message{Univariate{PointMass}}
isApplicable(::Type{SPClamp}, input_types::Vector{DataType}) = true