function matchPermutedCanonical(input_types::Vector{Type}, outbound_type::Type)
    # TODO: this implementation only works when the inbound types match the outbound type
    nothing_inputs = 0
    message_inputs = 0
    for input_type in input_types
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, outbound_type)
            message_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (message_inputs == 2)
end

mutable struct SPEqualityGaussian <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussian}) = Message{GaussianWeightedMeanPrecision}
isApplicable(::Type{SPEqualityGaussian}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Gaussian})

mutable struct SPEqualityGammaWishart <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGammaWishart}) = Message{Union{Gamma, Wishart}}
isApplicable(::Type{SPEqualityGammaWishart}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Union{Gamma, Wishart}})

mutable struct SPEqualityBernoulli <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityBernoulli}) = Message{Bernoulli}
isApplicable(::Type{SPEqualityBernoulli}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Bernoulli})

mutable struct SPEqualityBeta <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityBeta}) = Message{Beta}
isApplicable(::Type{SPEqualityBeta}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Beta})

mutable struct SPEqualityCategorical <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityCategorical}) = Message{Categorical}
isApplicable(::Type{SPEqualityCategorical}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Categorical})

mutable struct SPEqualityDirichlet <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityDirichlet}) = Message{Dirichlet}
isApplicable(::Type{SPEqualityDirichlet}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Dirichlet})

mutable struct SPEqualityPointMass <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityPointMass}) = Message{PointMass}
function isApplicable(::Type{SPEqualityPointMass}, input_types::Vector{Type})
    nothing_inputs = 0
    soft_inputs = 0
    point_mass_inputs = 0
    for input_type in input_types
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{SoftFactor})
            soft_inputs += 1
        elseif matches(input_type, Message{PointMass})
            point_mass_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (soft_inputs == 1) && (point_mass_inputs == 1)
end

mutable struct SPEqualityRGMP <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityRGMP}) = Message{Function}
isApplicable(::Type{SPEqualityRGMP}, input_types::Vector{Type}) = matchPermutedCanonical(input_types, Message{Function})

mutable struct SPEqualityGaussianRGMP <: SumProductRule{Equality} end
outboundType(::Type{SPEqualityGaussianRGMP}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPEqualityGaussianRGMP}, input_types::Vector{Type})
    nothing_inputs = 0
    function_inputs = 0
    gaussian_inputs = 0

    for input_type in input_types
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{Function})
            function_inputs += 1
        elseif matches(input_type, Message{Gaussian})
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (function_inputs == 1) && (gaussian_inputs == 1)
end
