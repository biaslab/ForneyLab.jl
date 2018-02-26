module NonlinearTest

using Base.Test
using ForneyLab
import ForneyLab: approximate, outboundType, isApplicable
import ForneyLab: SPNonlinearOutVGP, SPNonlinearIn1GVP, VBNonlinearOut, VBNonlinearIn1, VBNonlinearW, SVBNonlinearOutVGD, SVBNonlinearW, SVBNonlinearIn1GVD, MNonlinearGGD

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
g(x::Vector{Float64}) = x.^2 - 5.0
J_g(x::Vector{Float64}) = reshape(2*x, 1, 1)

@testset "SPNonlinearOutVGP" begin
    @test SPNonlinearOutVGP <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearOutVGP) == Message{Gaussian}
    @test isApplicable(SPNonlinearOutVGP, [Void, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPNonlinearOutVGP(nothing, Message(Univariate, Gaussian, m=2.0, v=3.0), Message(Univariate, PointMass, m=2.0), g, J_g) == Message(Univariate, Gaussian, m=-1.0, v=48.5)
    @test ruleSPNonlinearOutVGP(nothing, Message(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == Message(Multivariate, Gaussian, m=[-1.0], v=mat(48.5))
end

@testset "SPNonlinearIn1GVP" begin
    @test SPNonlinearIn1GVP <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearIn1GVP) == Message{Gaussian}
    @test isApplicable(SPNonlinearIn1GVP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPNonlinearIn1GVP(Message(Univariate, Gaussian, m=2.0, v=3.0), Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=2.0), g, J_g) == Message(Univariate, Gaussian, m=0.25*(2.0 + 9.0), v=0.25*(3.0 + 0.5)*0.25)
    @test ruleSPNonlinearIn1GVP(Message(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), Message(Multivariate, PointMass, m=[2.0]), Message(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == Message(Multivariate, Gaussian, m=[0.25*(2.0 + 9.0)], v=mat(0.25*(3.0 + 0.5)*0.25))
end

@testset "VBNonlinearOut" begin
    @test VBNonlinearOut <: NaiveVariationalRule{Nonlinear}
    @test outboundType(VBNonlinearOut) == Message{Gaussian}
    @test isApplicable(VBNonlinearOut, [Void, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBNonlinearOut(nothing, ProbabilityDistribution(Univariate, Gaussian, m=2.0, v=3.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), g, J_g) == Message(Univariate, Gaussian, m=-1.0, w=2.0)
    @test ruleVBNonlinearOut(nothing, ProbabilityDistribution(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == Message(Multivariate, Gaussian, m=[-1.0], w=mat(2.0))
end

@testset "VBNonlinearIn1" begin
    @test VBNonlinearIn1 <: NaiveVariationalRule{Nonlinear}
    @test outboundType(VBNonlinearIn1) == Message{Gaussian}
    @test isApplicable(VBNonlinearIn1, [ProbabilityDistribution, Void, ProbabilityDistribution]) 

    @test ruleVBNonlinearIn1(ProbabilityDistribution(Univariate, Gaussian, m=2.0, v=3.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), g, J_g) == Message(Univariate, Gaussian, m=11.0/4.0, w=32.0)
    @test ruleVBNonlinearIn1(ProbabilityDistribution(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), ProbabilityDistribution(Multivariate, PointMass, m=[2.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == Message(Multivariate, Gaussian, m=[11.0/4.0], w=mat(32.0))
end

@testset "VBNonlinearW" begin
    @test VBNonlinearW <: NaiveVariationalRule{Nonlinear}
    @test outboundType(VBNonlinearW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(VBNonlinearW, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBNonlinearW(ProbabilityDistribution(Univariate, Gaussian, m=4.0, v=5.0), ProbabilityDistribution(Univariate, Gaussian, m=2.0, v=3.0), nothing, g, J_g) == Message(Univariate, Gamma, a=1.5, b=0.5*(3.0*4.0^2 + 5.0 + (4.0*2.0 - 9.0 - 4.0)^2))
    @test ruleVBNonlinearW(ProbabilityDistribution(Multivariate, Gaussian, m=[4.0], v=mat(5.0)), ProbabilityDistribution(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), nothing, g, J_g) == Message(MatrixVariate, Wishart, v=mat(0.012820512820512818), nu=3.0)
end

@testset "SVBNonlinearOutVGD" begin
    @test SVBNonlinearOutVGD <: StructuredVariationalRule{Nonlinear}
    @test outboundType(SVBNonlinearOutVGD) == Message{Gaussian}
    @test isApplicable(SVBNonlinearOutVGD, [Void, Message{Gaussian}, ProbabilityDistribution]) 

    @test ruleSVBNonlinearOutVGD(nothing, Message(Univariate, Gaussian, m=2.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0), g, J_g) == Message(Univariate, Gaussian, m=-1.0, v=48.5)
    @test ruleSVBNonlinearOutVGD(nothing, Message(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0), g, J_g) == Message(Multivariate, Gaussian, m=[-1.0], v=mat(48.5))
end

@testset "SVBNonlinearIn1GVD" begin
    @test SVBNonlinearIn1GVD <: StructuredVariationalRule{Nonlinear}
    @test outboundType(SVBNonlinearIn1GVD) == Message{Gaussian}
    @test isApplicable(SVBNonlinearIn1GVD, [Message{Gaussian}, Void, ProbabilityDistribution]) 

    @test ruleSVBNonlinearIn1GVD(Message(Univariate, Gaussian, m=2.0, v=3.0), Message(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0), g, J_g) == Message(Univariate, Gaussian, m=0.25*(2.0 + 9.0), v=0.25*(3.0 + 0.5)*0.25)
    @test ruleSVBNonlinearIn1GVD(Message(Multivariate, Gaussian, m=[2.0], v=mat(3.0)), Message(Multivariate, PointMass, m=[2.0]), ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0), g, J_g) == Message(Multivariate, Gaussian, m=[0.25*(2.0 + 9.0)], v=mat(0.25*(3.0 + 0.5)*0.25))
end

@testset "SVBNonlinearW" begin
    @test SVBNonlinearW <: StructuredVariationalRule{Nonlinear}
    @test outboundType(SVBNonlinearW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(SVBNonlinearW, [ProbabilityDistribution, Void]) 

    @test ruleSVBNonlinearW(ProbabilityDistribution(Multivariate, Gaussian, m=[3.0, 2.0], v=[4.0 1.0; 1.0 3.0]), nothing, g, J_g) == Message(Univariate, Gamma, a=1.5, b=0.5*(4.0*3.0*4.0 - 2*4.0*1.0 + 4.0 + (4.0*2.0 - 9.0 - 3.0)^2))
    A = [2.0 2.0; 1.0 1.0] # Approximation for g_c at [2, 2]
    b = [-4.0, 0.0]
    @test ruleSVBNonlinearW(ProbabilityDistribution(Multivariate, Gaussian, m=[3.0, 4.0, 2.0, 2.0], v=[5.0 1.0 0.5 0.0; 1.0 4.0 2.0 0.5; 0.5 2.0 3.0 1.0; 0.0 0.5 1.0 2.0]), nothing, g_c, J_g_c) == Message(MatrixVariate, Wishart, v=cholinv( A*[3.0 1.0; 1.0 2.0]*A' - A*[0.5 2.0; 0.0 0.5] - [0.5 0.0; 2.0 0.5]*A' + [5.0 1.0; 1.0 4.0] + (A*[2.0, 2.0] + b - [3.0, 4.0])*(A*[2.0, 2.0] + b - [3.0, 4.0])'), nu=4.0)
end

@testset "MNonlinearGGD" begin
    @test MNonlinearGGD <: MarginalRule{Nonlinear}
    @test isApplicable(MNonlinearGGD, [Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution]) 

    @test ruleMNonlinearGGD(Message(Univariate, Gaussian, m=3.0, w=4.0), Message(Univariate, Gaussian, m=2.0, w=1.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0), g, J_g) == ProbabilityDistribution(Multivariate, Gaussian, m=[3.0, 2.0], w=[4.0 0.0; 0.0 1.0])*ProbabilityDistribution(Multivariate, Gaussian, m=[-9.0, 0.0], w=[2.0 -8.0; -8.0 32.0])
    @test ruleMNonlinearGGD(Message(Multivariate, Gaussian, m=[3.0], w=mat(4.0)), Message(Multivariate, Gaussian, m=[2.0], w=mat(1.0)), ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0), g, J_g) == ProbabilityDistribution(Multivariate, Gaussian, m=[3.0, 2.0], w=[4.0 0.0; 0.0 1.0])*ProbabilityDistribution(Multivariate, Gaussian, m=[-9.0, 0.0], w=[2.0 -8.0; -8.0 32.0])
end

@testset "averageEnergy" begin
    @test averageEnergy(Nonlinear, ProbabilityDistribution(Multivariate, Gaussian, m=[2.0], w=mat(2.0)), ProbabilityDistribution(Multivariate, PointMass, m=[2.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == 10.0723649429247
    @test averageEnergy(Nonlinear, ProbabilityDistribution(Univariate, Gaussian, m=2.0, w=2.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), g, J_g) == 10.0723649429247
    @test averageEnergy(Nonlinear, ProbabilityDistribution(Multivariate, Gaussian, m=[1.0, 2.0], v=[3.0 1.0; 1.0 2.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)), g, J_g) == 31.5723649429247
end

end # module