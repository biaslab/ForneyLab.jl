type SPConstant <: SumProductRule{Constant} end
outboundType(::Type{SPConstant}) = Message{PointMass}
isApplicable(::Type{SPConstant}, input_types) = true
# calculateMessage(::Type{SPConstant}, inputs) = ProbabilityDistribution{PointMass}(constant.value)