#------------------
# Sum-Product Rules
#------------------

mutable struct SPTransitionMixtureOutNCPX <: SumProductRule{TransitionMixture} end
outboundType(::Type{SPTransitionMixtureOutNCPX}) = Message{Categorical}
function isApplicable(::Type{SPTransitionMixtureOutNCPX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Nothing) || return false
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false
        else
            matches(input_type, Message{PointMass}) || return false
        end
    end
    return true
end

mutable struct SPTransitionMixtureIn1CNPX <: SumProductRule{TransitionMixture} end
outboundType(::Type{SPTransitionMixtureIn1CNPX}) = Message{Categorical}
function isApplicable(::Type{SPTransitionMixtureIn1CNPX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        elseif (i == 1)
            matches(input_type, Message{Categorical}) || return false
        else
            matches(input_type, Message{PointMass}) || return false
        end
    end
    return true
end


#------------
# Naive Rules
#------------

mutable struct VBTransitionMixtureZ <: NaiveVariationalRule{TransitionMixture} end
outboundType(::Type{VBTransitionMixtureZ}) = Message{Categorical}
function isApplicable(::Type{VBTransitionMixtureZ}, input_types::Vector{<:Type})
    (length(input_types) > 3) || return false
    for (i, input_type) in enumerate(input_types)
        if (i == 3)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct VBTransitionMixtureOut <: NaiveVariationalRule{TransitionMixture} end
outboundType(::Type{VBTransitionMixtureOut}) = Message{Categorical}
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

mutable struct VBTransitionMixtureIn1 <: NaiveVariationalRule{TransitionMixture} end
outboundType(::Type{VBTransitionMixtureIn1}) = Message{Categorical}
function isApplicable(::Type{VBTransitionMixtureIn1}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct VBTransitionMixtureA <: NaiveVariationalRule{TransitionMixture} end
outboundType(::Type{VBTransitionMixtureA}) = Message{Dirichlet}
function isApplicable(::Type{VBTransitionMixtureA}, input_types::Vector{<:Type})
    n_inputs = length(input_types)

    (Nothing_positions, p_positions) = matchPVInputs(input_types)
    n_Nothings = length(Nothing_positions)
    n_ps = length(p_positions)

    (n_Nothings == 1) || return false
    (n_Nothings + n_ps == n_inputs) || return false
    (1 in p_positions) || return false
    (2 in p_positions) || return false
    (3 in p_positions) || return false

    return true
end


#-----------------
# Structured Rules
#-----------------

mutable struct SVBTransitionMixtureZ <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureZ}) = Message{Categorical}
function isApplicable(::Type{SVBTransitionMixtureZ}, input_types::Vector{<:Type})
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

mutable struct SVBTransitionMixtureOutNCDX <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureOutNCDX}) = Message{Categorical}
function isApplicable(::Type{SVBTransitionMixtureOutNCDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Nothing) || return false
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false            
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct SVBTransitionMixtureIn1CNDX <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureIn1CNDX}) = Message{Categorical}
function isApplicable(::Type{SVBTransitionMixtureIn1CNDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        elseif (i == 1)
            matches(input_type, Message{Categorical}) || return false            
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct SVBTransitionMixtureA <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureA}) = Message{Dirichlet}
function isApplicable(::Type{SVBTransitionMixtureA}, input_types::Vector{<:Type})
    n_inputs = length(input_types)

    (Nothing_positions, p_positions) = matchPVInputs(input_types)
    n_Nothings = length(Nothing_positions)
    n_ps = length(p_positions)

    (n_Nothings == 1) || return false
    (n_Nothings + n_ps == n_inputs) || return false
    (1 in p_positions) || return false
    (2 in p_positions) || return false

    return true
end


#---------------
# Marginal Rules
#---------------

mutable struct MTransitionMixtureCCDX <: MarginalRule{TransitionMixture} end
function isApplicable(::Type{MTransitionMixtureCCDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false            
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end
