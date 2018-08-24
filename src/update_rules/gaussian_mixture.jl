mutable struct VBGaussianMixtureZBer <: NaiveVariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureZBer}) = Message{Bernoulli}
function isApplicable(::Type{VBGaussianMixtureZBer}, input_types::Vector{<:Type})
    (length(input_types) == 6) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct VBGaussianMixtureZCat <: NaiveVariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureZCat}) = Message{Categorical}
function isApplicable(::Type{VBGaussianMixtureZCat}, input_types::Vector{<:Type})
    (length(input_types) > 6) || return false
    iseven(length(input_types)) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

function matchPVInputs(input_types::Vector{<:Type})
    Nothing_positions = []
    p_positions = []
    for (i, input_type) in enumerate(input_types)
        if matches(input_type, ProbabilityDistribution)
            push!(p_positions, i)
        elseif (input_type == Nothing)
            push!(Nothing_positions, i)
        end
    end

    return (Nothing_positions, p_positions)
end

mutable struct VBGaussianMixtureM <: NaiveVariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureM}) = Message{GaussianMeanPrecision}
function isApplicable(::Type{VBGaussianMixtureM}, input_types::Vector{<:Type})
    n_inputs = length(input_types)
    iseven(n_inputs) || return false

    (Nothing_positions, p_positions) = matchPVInputs(input_types)
    n_Nothings = length(Nothing_positions)
    n_ps = length(p_positions)

    (n_Nothings == 1) || return false
    (n_Nothings + n_ps == n_inputs) || return false
    (1 in p_positions) || return false
    (2 in p_positions) || return false
    isodd(Nothing_positions[1]) || return false

    return true
end

mutable struct VBGaussianMixtureW <: NaiveVariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureW}) = Message{Union{Gamma, Wishart}}
function isApplicable(::Type{VBGaussianMixtureW}, input_types::Vector{<:Type})
    n_inputs = length(input_types)
    iseven(n_inputs) || return false

    (Nothing_positions, p_positions) = matchPVInputs(input_types)
    n_Nothings = length(Nothing_positions)
    n_ps = length(p_positions)

    (n_Nothings == 1) || return false
    (n_Nothings + n_ps == n_inputs) || return false
    (1 in p_positions) || return false
    (2 in p_positions) || return false
    iseven(Nothing_positions[1]) || return false

    return true
end

mutable struct VBGaussianMixtureOut <: NaiveVariationalRule{GaussianMixture} end
outboundType(::Type{VBGaussianMixtureOut}) = Message{GaussianWeightedMeanPrecision}
function isApplicable(::Type{VBGaussianMixtureOut}, input_types::Vector{<:Type})
    iseven(length(input_types)) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end