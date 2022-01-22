@sumProductRule(:node_type     => Delta{Extended},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPDeltaEOutNG)

@sumProductRule(:node_type     => Delta{Extended},
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPDeltaEIn1GG)

mutable struct SPDeltaEOutNGX <: SumProductRule{Delta{Extended}} end
outboundType(::Type{SPDeltaEOutNGX}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPDeltaEOutNGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end

mutable struct SPDeltaEInGX <: SumProductRule{Delta{Extended}} end
outboundType(::Type{SPDeltaEInGX}) = Message{Gaussian}
function isApplicable(::Type{SPDeltaEInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types
        if input_type == Nothing
            nothing_inputs += 1
        elseif input_type << Message{Gaussian}
            gaussian_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (gaussian_inputs == total_inputs - 1)
end

mutable struct MDeltaEInGX <: MarginalRule{Delta{Extended}} end
function isApplicable(::Type{MDeltaEInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
