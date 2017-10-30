type SPClamp <: SumProductRule{Clamp{Univariate{PointMass}}} end
outboundType(::Type{SPClamp}) = Message{Univariate{PointMass}}
isApplicable(::Type{SPClamp}, input_types::Vector{DataType}) = true