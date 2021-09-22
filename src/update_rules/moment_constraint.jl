mutable struct SPMomentConstraintOutG <: SumProductRule{MomentConstraint} end
outboundType(::Type{SPMomentConstraintOutG}) = Message{GaussianWeightedMeanPrecision}
isApplicable(::Type{SPMomentConstraintOutG}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
