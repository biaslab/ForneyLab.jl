@sumProductRule(:node_type     => Delta{Sampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message),
                :name          => SPDeltaSOutNM)

@sumProductRule(:node_type     => Delta{Sampling},
                :outbound_type => Message{Function},
                :inbound_types => (Message, Nothing),
                :name          => SPDeltaSIn1MN)

mutable struct SPDeltaSInGX <: SumProductRule{Delta{Sampling}} end
outboundType(::Type{SPDeltaSInGX}) = Message{Gaussian{Canonical}}
function isApplicable(::Type{SPDeltaSInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) && return false # Require any message on out

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

mutable struct SPDeltaSOutNMX <: SumProductRule{Delta{Sampling}} end
outboundType(::Type{SPDeltaSOutNMX}) = Message{SampleList}
function isApplicable(::Type{SPDeltaSOutNMX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false
    
    for input_type in input_types[2:end]
        (input_type << Message) || return false
    end

    return true
end

mutable struct SPDeltaSInMX <: SumProductRule{Delta{Sampling}} end
outboundType(::Type{SPDeltaSInMX}) = Message{Function}
function isApplicable(::Type{SPDeltaSInMX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) && return false

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

mutable struct MDeltaSInMGX <: MarginalRule{Delta{Sampling}} end
function isApplicable(::Type{MDeltaSInMGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
