abstract SPGaussianMeanVariance{PointMass, PointMass, Void} <: SumProductRule{GaussianMeanVariance}
outboundType(::Type{SPGaussianMeanVariance{PointMass, PointMass, Void}}) = ProbabilityDistribution{GaussianMeanVariance}
isApplicable(::Type{SPGaussianMeanVariance{PointMass, PointMass, Void}}, input_types) = (input_types[1] <: PointMass) && (input_types[2] <: PointMass)
# function calculateMessage(::Type{SPGaussianMeanVariance{PointMass, PointMass, Void}}, inputs)
#     m = unsafeMean(inputs[1])
#     v = unsafeMean(inputs[2])

#     return ProbabilityDistribution{GaussianMeanVariance}(m, v)
# end

abstract SPGaussianMeanVariance{Gaussian, PointMass, Void} <: SumProductRule{GaussianMeanVariance}
outboundType(::Type{SPGaussianMeanVariance{Gaussian, PointMass, Void}}) = ProbabilityDistribution{GaussianMeanVariance}
isApplicable(::Type{SPGaussianMeanVariance{Gaussian, PointMass, Void}}, input_types) = (input_types[1] <: Gaussian) && (input_types[2] <: PointMass)
# function calculateMessage(::Type{SPGaussianMeanVariance{Gaussian, PointMass, Void}}, inputs)
#     m = unsafeMean(inputs[1])
#     v = unsafeVar(inputs[1]) + unsafeMean(inputs[2])

#     return ProbabilityDistribution{GaussianMeanVariance}(m, v)
# end

