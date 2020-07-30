mutable struct SPPointMassConstraintOutG <: SumProductRule{PointMassConstraint} end
outboundType(::Type{SPPointMassConstraintOutG}) = Message{PointMass}
isApplicable(::Type{SPPointMassConstraintOutG}, input_types::Vector{<:Type}) = (length(input_types) == 1) && (input_types[1] == Nothing)
