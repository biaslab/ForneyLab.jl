#------------------
# Sum-Product Rules
#------------------

mutable struct SPTransitionMixtureOutNCCPX <: SumProductRule{TransitionMixture} end
outboundType(::Type{SPTransitionMixtureOutNCCPX}) = Message{Categorical}
function isApplicable(::Type{SPTransitionMixtureOutNCCPX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Nothing) || return false
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false
        elseif (i == 3)
            matches(input_type, Message{Categorical}) || return false
        else
            matches(input_type, Message{PointMass}) || return false
        end
    end
    return true
end

mutable struct SPTransitionMixtureIn1CNCPX <: SumProductRule{TransitionMixture} end
outboundType(::Type{SPTransitionMixtureIn1CNCPX}) = Message{Categorical}
function isApplicable(::Type{SPTransitionMixtureIn1CNCPX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        elseif (i == 1)
            matches(input_type, Message{Categorical}) || return false
        elseif (i == 3)
            matches(input_type, Message{Categorical}) || return false
        else
            matches(input_type, Message{PointMass}) || return false
        end
    end
    return true
end

mutable struct SPTransitionMixtureZCCNPX <: SumProductRule{TransitionMixture} end
outboundType(::Type{SPTransitionMixtureZCCNPX}) = Message{Categorical}
function isApplicable(::Type{SPTransitionMixtureZCCNPX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 3)
            (input_type == Nothing) || return false
        elseif (i == 1)
            matches(input_type, Message{Categorical}) || return false
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false
        else
            matches(input_type, Message{PointMass}) || return false
        end
    end
    return true
end


#-----------------
# Structured Rules
#-----------------

mutable struct SVBTransitionMixtureOutNCCDX <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureOutNCCDX}) = Message{Categorical}
function isApplicable(::Type{SVBTransitionMixtureOutNCCDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            (input_type == Nothing) || return false
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 3)
            matches(input_type, Message{Categorical}) || return false            
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct SVBTransitionMixtureIn1CNCDX <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureIn1CNCDX}) = Message{Categorical}
function isApplicable(::Type{SVBTransitionMixtureIn1CNCDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 2)
            (input_type == Nothing) || return false
        elseif (i == 1)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 3)
            matches(input_type, Message{Categorical}) || return false            
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct SVBTransitionMixtureZCCNDX <: StructuredVariationalRule{TransitionMixture} end
outboundType(::Type{SVBTransitionMixtureZCCNDX}) = Message{Categorical}
function isApplicable(::Type{SVBTransitionMixtureZCCNDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 3)
            (input_type == Nothing) || return false
        elseif (i == 1)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 2)
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

    return true
end


#---------------
# Marginal Rules
#---------------

mutable struct MTransitionMixtureCCCDX <: MarginalRule{TransitionMixture} end
function isApplicable(::Type{MTransitionMixtureCCCDX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 3)
            matches(input_type, Message{Categorical}) || return false            
        else
            matches(input_type, ProbabilityDistribution) || return false
        end
    end
    return true
end

mutable struct MTransitionMixtureCCCNX <: MarginalRule{TransitionMixture} end
function isApplicable(::Type{MTransitionMixtureCCCNX}, input_types::Vector{<:Type})
    for (i, input_type) in enumerate(input_types)
        if (i == 1)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 2)
            matches(input_type, Message{Categorical}) || return false            
        elseif (i == 3)
            matches(input_type, Message{Categorical}) || return false            
        else
            (input_type == Nothing) || return false
        end
    end
    return true
end
