mutable struct VBTransitionMixtureZ <: NaiveVariationalRule{TransitionMixture} end
outboundType(::Type{VBTransitionMixtureZ}) = Message{Categorical}
function isApplicable(::Type{VBTransitionMixtureZ}, input_types::Vector{<:Type})
    (length(input_types) > 2) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct VBTransitionMixtureOut <: NaiveVariationalRule{TransitionMixture} end
outboundType(::Type{VBTransitionMixtureOut}) = Message{Dirichlet}
function isApplicable(::Type{VBTransitionMixtureOut}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end