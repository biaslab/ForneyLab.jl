mutable struct SPCategoricalConstraintOut <: SumProductRule{CategoricalConstraint} end
outboundType(::Type{SPCategoricalConstraintOut}) = Message{Categorical}
isApplicable(::Type{SPCategoricalConstraintOut}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
