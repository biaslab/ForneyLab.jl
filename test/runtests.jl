# This file is automatically executed if you call Pkg.test("ForneyLab")

# include("testrunner.jl")

module ForneyLabTest

using Base.Test

@testset "ForneyLab" begin
    include("./test_helpers.jl")
    include("./test_dependency_graph.jl")
    include("./test_approximation.jl")
    include("./test_node.jl")
end

end