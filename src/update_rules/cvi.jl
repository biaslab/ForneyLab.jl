@sumProductRule(:node_type     => CVI,
                      :outbound_type => Message{SetSampleList},
                      :inbound_types => (Nothing, Message{FactorNode}),
                      :name          => SPCVIOutNFactorNode)

@sumProductRule(:node_type     => CVI,
                      :outbound_type => Message{FactorNode},
                      :inbound_types => (Message{FactorFunction}, Nothing),
                      :name          => SPCVIIn1Factor)

mutable struct SPCVIOutNFactorNodeX <: SumProductRule{CVI} end
outboundType(::Type{SPCVIOutNFactorNodeX}) = Message{SetSampleList}
function isApplicable(::Type{SPCVIOutNFactorNodeX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, Message{FactorNode}) || return false
    end

    return true
end

mutable struct SPCVIInFactorX <: SumProductRule{CVI} end
outboundType(::Type{SPCVIInFactorX}) = Message{FactorNode}
function isApplicable(::Type{SPCVIInFactorX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    nothing_inputs = 0
    factor_inputs = 0

    for input_type in input_types[1:end]
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, Message{FactorFunction})
            factor_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (factor_inputs == total_inputs-1)
end
