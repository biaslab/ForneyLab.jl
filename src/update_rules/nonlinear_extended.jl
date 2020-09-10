@sumProductRule(:node_type     => Nonlinear{Extended},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearEOutNG)

@sumProductRule(:node_type     => Nonlinear{Extended},
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPNonlinearEIn1GG)

mutable struct SPNonlinearEOutNGX <: SumProductRule{Nonlinear{Extended}} end
outboundType(::Type{SPNonlinearEOutNGX}) = Message{GaussianMeanVariance}
function isApplicable(::Type{SPNonlinearEOutNGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end

mutable struct SPNonlinearEInGX <: SumProductRule{Nonlinear{Extended}} end
outboundType(::Type{SPNonlinearEInGX}) = Message{Gaussian}
function isApplicable(::Type{SPNonlinearEInGX}, input_types::Vector{<:Type})
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

mutable struct MNonlinearEInGX <: MarginalRule{Nonlinear{Extended}} end
function isApplicable(::Type{MNonlinearEInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false # Indicates marginalization over outbound variable

    for input_type in input_types[2:end]
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end
