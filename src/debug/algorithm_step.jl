export AlgorithmStep, pushScore!, pushMessage!, pushMarginal!

struct AlgorithmStep
    messages::Vector{Vector{MessageData}}
    marginals::Vector{Vector{MarginalData}}
    score::Vector{ScoreData}
end

AlgorithmStep() = AlgorithmStep(Vector{Vector{MessageData}}(), Vector{Vector{MarginalData}}(), Vector{ScoreData}())

function pushScore!(step::AlgorithmStep, id::Union{String, Symbol}, value::Float64, type::String)
    push!(step.score, ScoreData(id, value, type))
    return step
end

function pushMessage!(step::AlgorithmStep, edgeID::Union{String, Symbol}, type::String, message)
    push!(step.messages[end], MessageData(edgeID, type, MessageSnapshot(message)))
    return step
end

function pushMarginal!(step::AlgorithmStep, id::Union{String, Symbol}, edgeIDs::Vector{String}, marginal)
    push!(step.marginals[end], MarginalData(id, edgeIDs, MarginalSnapshot(marginal)))
    return step
end
