module Julia

using ForneyLab

export
messagePassingAlgorithm,
rule

include("message_passing.jl")
include("update_rules/gaussian_mean_variance.jl")

end