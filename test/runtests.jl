# This file is automatically executed if you call Pkg.test("ForneyLab")

module ForneyLabTest

using Test

@testset "ForneyLab" begin
    include("./test_helpers.jl")
    include("./test_factor_node.jl")
    include("./test_interface.jl")
    include("./test_probability_distribution.jl")
    include("./test_edge.jl")
    include("./test_variable.jl")
    include("./test_factor_graph.jl")

    # Factor nodes
    include("factor_nodes/test_clamp.jl")
    include("factor_nodes/test_equality.jl")
    include("factor_nodes/test_addition.jl")
    include("factor_nodes/test_multiplication.jl")
    include("factor_nodes/test_exponential.jl")
    include("factor_nodes/test_gaussian_mean_variance.jl")
    include("factor_nodes/test_gaussian_mean_precision.jl")
    include("factor_nodes/test_gaussian_weighted_mean_precision.jl")
    include("factor_nodes/test_gaussian.jl")
    include("factor_nodes/test_gamma.jl")
    include("factor_nodes/test_log_normal.jl")
    include("factor_nodes/test_wishart.jl")
    include("factor_nodes/test_bernoulli.jl")
    include("factor_nodes/test_categorical.jl")
    include("factor_nodes/test_contingency.jl")
    include("factor_nodes/test_transition.jl")
    include("factor_nodes/test_beta.jl")
    include("factor_nodes/test_dirichlet.jl")
    include("factor_nodes/test_gaussian_mixture.jl")
    include("factor_nodes/test_probit.jl")
    include("factor_nodes/test_logit.jl")
    include("factor_nodes/test_softmax.jl")
    include("factor_nodes/test_dot_product.jl")
    include("factor_nodes/test_poisson.jl")
    include("factor_nodes/test_sample_list.jl")
    include("factor_nodes/test_nonlinear_unscented.jl")
    include("factor_nodes/test_nonlinear_sampling.jl")
    include("factor_nodes/test_nonlinear_extended.jl")
    include("./test_dependency_graph.jl")
    include("./test_message_passing.jl")
    include("./test_marginals.jl")

    # Algorithms
    include("./algorithms/test_cluster.jl")
    include("./algorithms/test_posterior_factor.jl")
    include("./algorithms/test_posterior_factorization.jl")
    include("./algorithms/test_inference_algorithm.jl")
    include("./algorithms/test_joint_marginals.jl")
    include("./algorithms/test_sum_product.jl")
    include("./algorithms/test_naive_variational_bayes.jl")
    include("./algorithms/test_structured_variational_bayes.jl")
    include("./algorithms/test_expectation_propagation.jl")

    # Assemblers
    include("engines/test_message_passing_assemblers.jl")
    include("engines/test_free_energy_assemblers.jl")

    # Engines
    include("./engines/julia/test_generators.jl")

    # Composite node
    include("factor_nodes/test_composite.jl")
end

end
