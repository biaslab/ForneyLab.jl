mutable struct SPPointMassConstraintOut <: SumProductRule{PointMassConstraint} end
outboundType(::Type{SPPointMassConstraintOut}) = Message{PointMass}
isApplicable(::Type{SPPointMassConstraintOut}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
