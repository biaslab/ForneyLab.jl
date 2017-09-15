type SPEqualityGaussianMV <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussianMV}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPEqualityGaussianMV}, input_types::Vector{DataType})
    sum(input_types .== Message{GaussianMeanVariance}) >= length(input_types) - 1
end

type SPEqualityGamma <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGamma}) = Message{Gamma}
function isApplicable(::Type{SPEqualityGamma}, input_types::Vector{DataType})
    sum(input_types .== Message{Gamma}) >= length(input_types) - 1
end