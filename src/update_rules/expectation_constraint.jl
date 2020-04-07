mutable struct SPExpectationConstraintOutG <: SumProductRule{ExpectationConstraint} end
outboundType(::Type{SPExpectationConstraintOutG}) = Message{GaussianWeightedMeanPrecision}
isApplicable(::Type{SPExpectationConstraintOutG}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
