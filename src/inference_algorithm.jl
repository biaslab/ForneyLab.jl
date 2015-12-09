export InferenceAlgorithm, GenericAlgorithm

abstract InferenceAlgorithm

type GenericAlgorithm <: InferenceAlgorithm
    execute::Function
end

Base.deepcopy(::InferenceAlgorithm) = error("deepcopy(::InferenceAlgorithm) is not possible. You should construct a new InferenceAlgorithm.")