module ApproximationTest

using Base.Test
import ForneyLab: Approximation, Laplace, GaussianMeanVariance

@testset "Approximation" begin
    # parameters of approximations should be subtype-constrained
    @test Approximation{GaussianMeanVariance, Laplace} <: Approximation
    @test_throws Exception Approximation{GaussianMeanVariance, Int64}
    @test_throws Exception Approximation{Int64, Laplace}
end

end #module