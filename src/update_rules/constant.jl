type SPConstant <: SumProductRule{Constant} end
outboundType(::Type{SPConstant}) = Message{PointMass}
isApplicable(::Type{SPConstant}, input_types::Vector{DataType}) = true