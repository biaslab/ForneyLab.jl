@sumProductRule(:node_type     => Delta{Conjugate},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message),
                :name          => SPDeltaCOutNM)

@sumProductRule(:node_type     => Delta{Conjugate},
                :outbound_type => Message{FactorNode},
                :inbound_types => (Message, Nothing),
                :name          => SPDeltaCIn1MN)

mutable struct SPDeltaCInGX <: SumProductRule{Delta{Conjugate}} end
outboundType(::Type{SPDeltaCInGX}) = Message{GaussianWeightedMeanPrecision}
function isApplicable(::Type{SPDeltaCInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false # Require any message on out

    nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types[2:end]
        if input_type == Nothing
            nothing_inputs += 1
        elseif input_type << Message{Gaussian}
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs == total_inputs - 2)
end

mutable struct SPDeltaCOutNMX <: SumProductRule{Delta{Conjugate}} end
outboundType(::Type{SPDeltaCOutNMX}) = Message{SampleList}
function isApplicable(::Type{SPDeltaCOutNMX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false
    
    for input_type in input_types[2:end]
        (input_type << Message) || return false
    end

    return true
end

mutable struct SPDeltaCInMX <: SumProductRule{Delta{Conjugate}} end
outboundType(::Type{SPDeltaCInMX}) = Message{FactorNode}
function isApplicable(::Type{SPDeltaCInMX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types[2:end]
        if input_type == Nothing
            nothing_inputs += 1
        elseif input_type << Message{Gaussian}
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs != total_inputs - 2) # Rule does not apply if all inbounds are Gaussian
end

mutable struct MDeltaCInMGX <: MarginalRule{Delta{Conjugate}} end
function isApplicable(::Type{MDeltaCInMGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
