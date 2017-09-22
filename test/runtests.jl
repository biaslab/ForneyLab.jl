# This file is automatically executed if you call Pkg.test("ForneyLab")

# include("testrunner.jl")

module ForneyLabTest

using Base.Test

@testset "ForneyLab" begin
    include("./test_helpers.jl")
    include("./test_factor_node.jl")
    include("./test_interface.jl")
    include("./test_probability_distribution.jl")
    include("./test_approximation.jl")
    include("./test_edge.jl")
    include("./test_variable.jl")
    include("./test_factor_graph.jl")

    # Factor nodes
    include("factor_nodes/test_clamp.jl")
    include("factor_nodes/test_equality.jl")
    include("factor_nodes/test_addition.jl")
    include("factor_nodes/test_multiplication.jl")
    include("factor_nodes/test_gaussian.jl")
    include("factor_nodes/test_gaussian_mean_variance.jl")
    include("factor_nodes/test_gaussian_mean_precision.jl")
    include("factor_nodes/test_gamma.jl")

    include("./test_dependency_graph.jl")
    include("./test_message_passing.jl")

    # Algorithms
    include("./algorithms/sum_product/test_sum_product.jl")

    # Engines
    include("./engines/julia/test_message_passing.jl")

    # Composite node
    include("factor_nodes/test_composite.jl")
end

end