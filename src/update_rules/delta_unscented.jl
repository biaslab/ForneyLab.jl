@sumProductRule(:node_type     => Delta{Unscented},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPDeltaUTOutNG)

@sumProductRule(:node_type     => Delta{Unscented},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPDeltaUTIn1GG)

mutable struct SPDeltaUTOutNGX <: SumProductRule{Delta{Unscented}} end
outboundType(::Type{SPDeltaUTOutNGX}) = Message{Gaussian{Moments}}
function isApplicable(::Type{SPDeltaUTOutNGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end

mutable struct SPDeltaUTInGX <: SumProductRule{Delta{Unscented}} end
outboundType(::Type{SPDeltaUTInGX}) = Message{Gaussian}
function isApplicable(::Type{SPDeltaUTInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{Gaussian})
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs == total_inputs - 1)
end

mutable struct MDeltaUTInGX <: MarginalRule{Delta{Unscented}} end
function isApplicable(::Type{MDeltaUTInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
