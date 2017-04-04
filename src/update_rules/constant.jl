abstract SPConstant <: SumProductRule{Constant}
outboundType(::Type{SPConstant}) = ProbabilityDistribution{PointMass}
isApplicable(::Type{SPConstant}, input_types) = true
# calculateMessage(::Type{SPConstant}, inputs) = ProbabilityDistribution{PointMass}(constant.value)