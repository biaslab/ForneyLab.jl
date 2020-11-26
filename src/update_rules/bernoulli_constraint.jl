mutable struct SPBernoulliConstraintOut <: SumProductRule{BernoulliConstraint} end
outboundType(::Type{SPBernoulliConstraintOut}) = Message{Bernoulli}
isApplicable(::Type{SPBernoulliConstraintOut}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
