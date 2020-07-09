export AlgorithmStep

mutable struct AlgorithmStep
    messages::Vector{Vector{MessageData}}
    marginals::Vector{Vector{MarginalData}}
    score::Vector{ScoreData}
end
