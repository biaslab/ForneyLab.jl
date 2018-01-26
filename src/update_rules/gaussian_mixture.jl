mutable struct VBGaussianMixtureZBer <: VariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureZBer}) = Message{Bernoulli}
function isApplicable(::Type{VBGaussianMixtureZBer}, input_types::Vector{<:Type})
    (length(input_types) == 6) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Void) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct VBGaussianMixtureZCat <: VariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureZCat}) = Message{Categorical}
function isApplicable(::Type{VBGaussianMixtureZCat}, input_types::Vector{<:Type})
    (length(input_types) > 6) || return false
    iseven(length(input_types)) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Void) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

function matchPVInputs(input_types::Vector{<:Type})
    void_positions = []
    p_positions = []
    for (i, input_type) in enumerate(input_types)
        if matches(input_type, ProbabilityDistribution)
            push!(p_positions, i)
        elseif (input_type == Void)
            push!(void_positions, i)
        end
    end

    return (void_positions, p_positions)
end

mutable struct VBGaussianMixtureM <: VariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureM}) = Message{Gaussian}
function isApplicable(::Type{VBGaussianMixtureM}, input_types::Vector{<:Type})
    n_inputs = length(input_types)
    iseven(n_inputs) || return false
    
    (void_positions, p_positions) = matchPVInputs(input_types)
    n_voids = length(void_positions)
    n_ps = length(p_positions)

    (n_voids == 1) || return false
    (n_voids + n_ps == n_inputs) || return false
    (1 in p_positions) || return false
    (2 in p_positions) || return false
    isodd(void_positions[1]) || return false
    
    return true
end

mutable struct VBGaussianMixtureW <: VariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureW}) = Message{Union{Gamma, Wishart}}
function isApplicable(::Type{VBGaussianMixtureW}, input_types::Vector{<:Type})
    n_inputs = length(input_types)
    iseven(n_inputs) || return false
    
    (void_positions, p_positions) = matchPVInputs(input_types)
    n_voids = length(void_positions)
    n_ps = length(p_positions)

    (n_voids == 1) || return false
    (n_voids + n_ps == n_inputs) || return false
    (1 in p_positions) || return false
    (2 in p_positions) || return false
    iseven(void_positions[1]) || return false
    
    return true
end

mutable struct VBGaussianMixtureOut <: VariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureOut}) = Message{Gaussian}
function isApplicable(::Type{VBGaussianMixtureOut}, input_types::Vector{<:Type})
    iseven(length(input_types)) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Void) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end