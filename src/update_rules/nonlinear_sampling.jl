@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message),
                :name          => SPNonlinearSOutNM)

@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{Function},
                :inbound_types => (Message, Nothing),
                :name          => SPNonlinearSIn1MN)

mutable struct SPNonlinearSOutNGX <: SumProductRule{Nonlinear{Sampling}} end
outboundType(::Type{SPNonlinearSOutNGX}) = Message{SampleList}
function isApplicable(::Type{SPNonlinearSOutNGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end

mutable struct SPNonlinearSInGX <: SumProductRule{Nonlinear{Sampling}} end
outboundType(::Type{SPNonlinearSInGX}) = Message{GaussianWeightedMeanPrecision}
function isApplicable(::Type{SPNonlinearSInGX}, input_types::Vector{<:Type})
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

mutable struct MNonlinearSInGX <: MarginalRule{Nonlinear{Sampling}} end
function isApplicable(::Type{MNonlinearSInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
