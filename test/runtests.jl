# This file is automatically executed if you call Pkg.test("ForneyLab")

# include("testrunner.jl")

module ForneyLabTest

using Base.Test

@testset "ForneyLab" begin
    include("./test_helpers.jl")
    include("./test_factor_node.jl")
    include("./test_interface.jl")
    # include("./test_message.jl")
    include("./test_probability_distribution.jl")
    include("./test_approximation.jl")
    include("./test_edge.jl")
    include("./test_variable.jl")
    include("./test_factor_graph.jl")

    # Factor nodes
    # include("factor_nodes/constant.jl")
    # include("factor_nodes/equality.jl")
    # include("factor_nodes/gaussian_mean_variance.jl")

    # include("./test_dependency_graph.jl")
end

end