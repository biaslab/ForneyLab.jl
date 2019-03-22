module NonlinearTest

using Test
using ForneyLab
import ForneyLab: approximate, outboundType, isApplicable
import ForneyLab: SPNonlinearOutVG, SPNonlinearIn1GV

g_a(x) = x.^2 - 5.0
g_prime_a(x) = 2*x

g_b(x) = 2.0*x[1].^2 + 2.0*x[2].^2 + x[1]*x[2] - 5.0
J_g_b(x) = [4.0*x[1] + x[2] 4.0*x[2] + x[1]]

g_c(x) = [x[1]*x[2], x[1] + x[2]]
J_g_c(x) = [x[2] x[1]; 1.0 1.0]

@testset "approximate" begin
    # g: R -> R
    x_hat = 2.0
    @test ForneyLab.approximate(x_hat, g_a, g_prime_a) == (4.0, -9.0)

    # g: R^2 -> R
    x_hat = [2.0, 2.0]
    @test ForneyLab.approximate(x_hat, g_b, J_g_b) == ([10.0 10.0], [-25.0])

    # g: R^2 -> R^2
    x_hat = [2.0, 2.0]
    @test ForneyLab.approximate(x_hat, g_c, J_g_c) == ([2.0 2.0; 1.0 1.0], [-4.0, 0.0])
end

#-------------
# Update rules
#-------------

g(x::Float64) = x^2 - 5.0
J_g(x::Float64) = 2*x
g(x::Vector{Float64}) = x.^2 .- 5.0
J_g(x::Vector{Float64}) = reshape(2*x, 1, 1)

@testset "SPNonlinearOutVG" begin
    @test SPNonlinearOutVG <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearOutVG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearOutVG, [Nothing, Message{Gaussian}])

    @test ruleSPNonlinearOutVG(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), g, J_g) == Message(Univariate, GaussianMeanVariance, m=-1.0, v=48.0)
    @test ruleSPNonlinearOutVG(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), g, J_g) == Message(Multivariate, GaussianMeanVariance, m=[-1.0], v=mat(48.0))
end

@testset "SPNonlinearIn1GV" begin
    @test SPNonlinearIn1GV <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearIn1GV) == Message{GaussianMeanPrecision}
    @test isApplicable(SPNonlinearIn1GV, [Message{Gaussian}, Nothing])

    @test ruleSPNonlinearIn1GV(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, PointMass, m=2.0), g, J_g) == Message(Univariate, GaussianMeanPrecision, m=0.25*(2.0 + 9.0), w=1.0/(0.25*3.0*0.25))
    @test ruleSPNonlinearIn1GV(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, PointMass, m=[2.0]), g, J_g) == Message(Multivariate, GaussianMeanPrecision, m=[0.25*(2.0 + 9.0)], w=mat(1.0/(0.25*3.0*0.25)))
end

end # module
