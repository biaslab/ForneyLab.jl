mutable struct SPGaussianMeanVarianceConstraintOut <: SumProductRule{GaussianMeanVarianceConstraint} end
outboundType(::Type{SPGaussianMeanVarianceConstraintOut}) = Message{GaussianWeightedMeanPrecision}
isApplicable(::Type{SPGaussianMeanVarianceConstraintOut}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
