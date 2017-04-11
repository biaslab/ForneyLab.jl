type SPGaussianMeanVariancePPV <: SumProductRule{GaussianMeanVariance} end
outboundType(::Type{SPGaussianMeanVariancePPV}) = Message{GaussianMeanVariance}
isApplicable(::Type{SPGaussianMeanVariancePPV}, input_types::Vector{DataType}) = (input_types[1] == Message{PointMass}) && (input_types[2] == Message{PointMass})

# TODO: this function will not work for a GaussianMeanPrecision input on interface 1
type SPGaussianMeanVarianceGPV <: SumProductRule{GaussianMeanVariance} end
outboundType(::Type{SPGaussianMeanVarianceGPV}) = Message{GaussianMeanVariance}
isApplicable(::Type{SPGaussianMeanVarianceGPV}, input_types::Vector{DataType}) = (input_types[1] == Message{GaussianMeanVariance}) && (input_types[2] == Message{PointMass})
