# Julia-specific generators
include("generators.jl")

# Julia-specific update implementations
include("update_rules/equality.jl")
include("update_rules/addition.jl")
include("update_rules/multiplication.jl")
include("update_rules/exponential.jl")
include("update_rules/gaussian_mean_variance.jl")
include("update_rules/gaussian_mean_precision.jl")
include("update_rules/gaussian_weighted_mean_precision.jl")
include("update_rules/gamma.jl")
include("update_rules/log_normal.jl")
include("update_rules/wishart.jl")
include("update_rules/bernoulli.jl")
include("update_rules/categorical.jl")
include("update_rules/transition.jl")
include("update_rules/beta.jl")
include("update_rules/dirichlet.jl")
include("update_rules/gaussian_mixture.jl")
include("update_rules/transition_mixture.jl")
include("update_rules/probit.jl")
include("update_rules/logit.jl")
include("update_rules/softmax.jl")
include("update_rules/nonlinear_unscented.jl")
include("update_rules/nonlinear_sampling.jl")
include("update_rules/nonlinear_extended.jl")
include("update_rules/nonlinear_conjugate.jl")
include("update_rules/dot_product.jl")
include("update_rules/poisson.jl")
include("update_rules/moment_constraint.jl")
include("update_rules/chance_constraint.jl")
include("update_rules/point_mass_constraint.jl")
include("update_rules/sample_list.jl")