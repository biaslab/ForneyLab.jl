mutable struct SPChanceConstraintOutG <: SumProductRule{ChanceConstraint} end
outboundType(::Type{SPChanceConstraintOutG}) = Message{Gaussian{Canonical}}
isApplicable(::Type{SPChanceConstraintOutG}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
