module ApproximationTest

using Base.Test
import ForneyLab: Approximation, Laplace, Gaussian

@testset "Approximation" begin
    # parameters of approximations should be subtype-constrained
    @test Approximation{Gaussian, Laplace} == Approximation{Gaussian, Laplace}
    @test_throws Exception Approximation{Gaussian, Int64}
    @test_throws Exception Approximation{Int64, Laplace}
end

end #module