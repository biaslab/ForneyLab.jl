abstract SPGaussianMeanVariancePPV <: SumProductRule{GaussianMeanVariance}
outboundType(::Type{SPGaussianMeanVariancePPV}) = ProbabilityDistribution{GaussianMeanVariance}
isApplicable(::Type{SPGaussianMeanVariancePPV}, input_types) = (input_types[1] <: PointMass) && (input_types[2] <: PointMass)
# function calculateMessage(::Type{SPGaussianMeanVariancePPV}, inputs)
#     m = unsafeMean(inputs[1])
#     v = unsafeMean(inputs[2])

#     return ProbabilityDistribution{GaussianMeanVariance}(m, v)
# end

abstract SPGaussianMeanVarianceGPV <: SumProductRule{GaussianMeanVariance}
outboundType(::Type{SPGaussianMeanVarianceGPV}) = ProbabilityDistribution{GaussianMeanVariance}
isApplicable(::Type{SPGaussianMeanVarianceGPV}, input_types) = (input_types[1] <: Gaussian) && (input_types[2] <: PointMass)
# function calculateMessage(::Type{SPGaussianMeanVarianceGPV}, inputs)
#     m = unsafeMean(inputs[1])
#     v = unsafeVar(inputs[1]) + unsafeMean(inputs[2])

#     return ProbabilityDistribution{GaussianMeanVariance}(m, v)
# end

