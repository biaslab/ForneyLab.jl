@sumProductRule(:node_type     => Nonlinear{Unscented},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearUTOutNG)

@sumProductRule(:node_type     => Nonlinear{Unscented},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPNonlinearUTIn1GG)

mutable struct SPNonlinearUTOutNGX <: SumProductRule{Nonlinear{Unscented}} end
outboundType(::Type{SPNonlinearUTOutNGX}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPNonlinearUTOutNGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end

mutable struct SPNonlinearUTInGX <: SumProductRule{Nonlinear{Unscented}} end
outboundType(::Type{SPNonlinearUTInGX}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPNonlinearUTInGX}, input_types::Vector{<:Type})
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

    return (nothing_inputs == 1) && (gaussian_inputs == total_inputs-1)
end

mutable struct MNonlinearUTInGX <: MarginalRule{Nonlinear{Unscented}} end
function isApplicable(::Type{MNonlinearUTInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
