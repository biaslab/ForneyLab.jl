@sumProductRule(:node_type     => Nonlinear{Unscented},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearUTOutNG)

@sumProductRule(:node_type     => Nonlinear{Unscented},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPNonlinearUTIn1GG)

mutable struct ruleSPNonlinearUTOutNGX <: SumProductRule{Nonlinear{Unscented}} end
outboundType(::Type{ruleSPNonlinearUTOutNGX}) = Message{GaussianMeanVariance}
function isApplicable(::Type{ruleSPNonlinearUTOutNGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types
        matches(input_type, Message{Gaussian}) || return false
    end

    return true
end

mutable struct ruleSPNonlinearUTInGX <: SumProductRule{Nonlinear{Unscented}} end
outboundType(::Type{ruleSPNonlinearUTInGX}) = Message{GaussianMeanVariance}
function isApplicable(::Type{ruleSPNonlinearUTInGX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    Nothing_inputs = 0
    gaussian_inputs = 0
    for input_type in input_types
        if input_type == Nothing
            Nothing_inputs += 1
        elseif matches(input_type, Message{Gaussian})
            gaussian_inputs += 1
        end
    end

    return (Nothing_inputs == 1) && (gaussian_inputs == total_inputs-1)
end

@sumProductRule(:node_type     => Nonlinear{ImportanceSampling},
                :outbound_type => Message{Function},
                :inbound_types => (Message{FactorFunction}, Nothing),
                :name          => SPNonlinearISIn1MN)

@sumProductRule(:node_type     => Nonlinear{ImportanceSampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearISOutNG)