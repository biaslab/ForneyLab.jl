export messagePassingAlgorithm

include("message_passing.jl")
include("sum_product.jl")
include("variational_bayes.jl")
include("expectation_propagation.jl")

include("update_rules/equality.jl")
include("update_rules/addition.jl")
include("update_rules/multiplication.jl")
include("update_rules/exponential.jl")
include("update_rules/gaussian_mean_variance.jl")
include("update_rules/gaussian_mean_precision.jl")
include("update_rules/gamma.jl")
include("update_rules/log_normal.jl")
include("update_rules/wishart.jl")
include("update_rules/bernoulli.jl")
include("update_rules/categorical.jl")
include("update_rules/transition.jl")
include("update_rules/beta.jl")
include("update_rules/dirichlet.jl")
include("update_rules/gaussian_mixture.jl")
include("update_rules/sigmoid.jl")
include("update_rules/nonlinear.jl")
include("update_rules/dot_product.jl")
include("update_rules/gain_equality.jl")
