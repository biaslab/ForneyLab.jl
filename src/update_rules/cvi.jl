@naiveVariationalRule(:node_type     => CVI,
                      :outbound_type => Message{SetSampleList},
                      :inbound_types => (Nothing, ProbabilityDistribution),
                      :name          => VBCVIOutVD)

@naiveVariationalRule(:node_type     => CVI,
                      :outbound_type => Message{FactorNode},
                      :inbound_types => (ProbabilityDistribution, Nothing),
                      :name          => VBCVIIn1MV)

mutable struct VBCVIOutVDX <: NaiveVariationalRule{CVI} end
outboundType(::Type{VBCVIOutVDX}) = Message{SetSampleList}
function isApplicable(::Type{VBCVIOutVDX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] == Nothing) || return false

    for input_type in input_types[2:end]
        matches(input_type, ProbabilityDistribution) || return false
    end

    return true
end

mutable struct VBCVIInX <: NaiveVariationalRule{CVI} end
outboundType(::Type{VBCVIInX}) = Message{FactorNode}
function isApplicable(::Type{VBCVIInX}, input_types::Vector{<:Type})
    total_inputs = length(input_types)
    (total_inputs > 2) || return false
    (input_types[1] != Nothing) || return false

    nothing_inputs = 0
    prob_inputs = 0

    for input_type in input_types[1:end]
        if input_type == Nothing
            nothing_inputs += 1
        elseif matches(input_type, ProbabilityDistribution)
            prob_inputs += 1
        end
    end

    return (nothing_inputs == 1) && (prob_inputs == total_inputs-1)
end
