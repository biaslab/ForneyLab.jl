export AlgorithmStep

mutable struct AlgorithmStep
    messages::Vector{Vector{MessageData}}
    marginals::Vector{Vector{MarginalData}}
    score::Vector{ScoreData}
end

function AlgorithmStep()
    step = AlgorithmStep(Vector{Vector{MessageData}}(), Vector{Vector{MarginalData}}(), Vector{ScoreData}())
    push!(step.messages, Vector{MessageData}())
    push!(step.marginals, Vector{MarginalData}())
    return step
end
