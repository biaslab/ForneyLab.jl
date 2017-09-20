type SPEqualityGaussian <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussian}) = Message{Gaussian}
function isApplicable(::Type{SPEqualityGaussian}, input_types::Vector{DataType})
    void_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types
        if input_type == Void
            void_inputs += 1
        elseif input_type == Message{Gaussian}
            gaussian_inputs += 1
        end
    end
    
    return (void_inputs == 1) && (gaussian_inputs == 2)
end