module NonlinearTest

using Base.Test
using ForneyLab
import ForneyLab: approximate, outboundType, isApplicable
import ForneyLab: SPNonlinearOutVGP

@testset "approximate" begin
    # g: R -> R
    x_hat = 2.0
    g_a(x) = x.^2 - 5.0
    g_prime_a(x) = 2*x
    @test ForneyLab.approximate(x_hat, g_a, g_prime_a) == (4.0, -9.0)

    # g: R^2 -> R
    x_hat = [2.0, 2.0]
    g_b(x) = 2.0*x[1].^2 + 2.0*x[2].^2 + x[1]*x[2] - 5.0
    J_g_b(x) = [4.0*x[1] + x[2] 4.0*x[2] + x[1]]
    @test ForneyLab.approximate(x_hat, g_b, J_g_b) == ([10.0 10.0], [-25.0])

    # g: R^2 -> R^2
    x_hat = [2.0, 2.0]
    g_c(x) = [x[1]*x[2], x[1] + x[2]]
    J_g_c(x) = [x[2] x[1]; 1.0 1.0]
    @test ForneyLab.approximate(x_hat, g_c, J_g_c) == ([2.0 2.0; 1.0 1.0], [-4.0, 0.0])
end

#-------------
# Update rules
#-------------

g(x::Float64) = x^2 - 5.0
J_g(x::Float64) = 2*x
g(x::Vector{Float64}) = x.^2 - 5.0
J_g(x::Vector{Float64}) = reshape(2*x, 1, 1)

@testset "SPNonlinearOutVGP" begin
    @test SPNonlinearOutVGP <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearOutVGP) == Message{Gaussian}
    @test isApplicable(SPNonlinearOutVGP, [Void, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPNonlinearOutVGP(nothing, Message(Univariate, Gaussian, m=2.0, v=3.0), Message(Univariate, PointMass, m=2.0), g, J_g) == Message(Univariate, Gaussian, m=-1.0, v=48.5)
    @test ruleSPNonlinearOutVGP(nothing, Message(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == Message(Multivariate, Gaussian, m=[-1.0], v=mat(48.5))
end

end # module