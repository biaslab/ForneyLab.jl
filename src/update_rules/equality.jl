type SPEqualityGaussian <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussian}) = Message{GaussianWeightedMeanPrecision}
function isApplicable(::Type{SPEqualityGaussian}, input_types::Vector{DataType})
    void_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types
        if input_type == Void
            void_inputs += 1
        elseif nodeType(input_type) <: Gaussian
            gaussian_inputs += 1
        end
    end
    
    return (void_inputs == 1) && (gaussian_inputs == 2)
end

nodeType{T}(::Type{Message{T}}) = T