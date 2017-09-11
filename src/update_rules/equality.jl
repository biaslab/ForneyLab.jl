type SPEqualityGaussianMV <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussianMV}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPEqualityGaussianMV}, input_types::Vector{DataType})
    sum(input_types .== Message{GaussianMeanVariance}) >= length(input_types) - 1
end