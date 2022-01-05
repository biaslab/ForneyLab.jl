@sumProductRule(:node_type     => Nonlinear{Conjugate},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message),
                :name          => SPNonlinearCOutNM)

@sumProductRule(:node_type     => Nonlinear{Conjugate},
                :outbound_type => Message{FactorNode},
                :inbound_types => (Message, Nothing),
                :name          => SPNonlinearCIn1MN)

mutable struct SPNonlinearCInGX <: SumProductRule{Nonlinear{Conjugate}} end
outboundType(::Type{SPNonlinearCInGX}) = Message{GaussianWeightedMeanPrecision}
function isApplicable(::Type{SPNonlinearCInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false # Require any message on out

    nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types[2:end]
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{Gaussian})
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs == total_inputs - 2)
end

mutable struct SPNonlinearCOutNMX <: SumProductRule{Nonlinear{Conjugate}} end
outboundType(::Type{SPNonlinearCOutNMX}) = Message{SampleList}
function isApplicable(::Type{SPNonlinearCOutNMX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false
    
    for input_type in input_types[2:end]
        matches(input_type, Message) || return false
    end

    return true
end

mutable struct SPNonlinearCInMX <: SumProductRule{Nonlinear{Conjugate}} end
outboundType(::Type{SPNonlinearCInMX}) = Message{FactorNode}
function isApplicable(::Type{SPNonlinearCInMX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types[2:end]
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{Gaussian})
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs != total_inputs - 2) # Rule does not apply if all inbounds are Gaussian
end

mutable struct MNonlinearCInMGX <: MarginalRule{Nonlinear{Conjugate}} end
function isApplicable(::Type{MNonlinearCInMGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
