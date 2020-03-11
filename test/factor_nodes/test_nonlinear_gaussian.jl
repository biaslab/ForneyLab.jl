module NonlinearGaussianTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights
using ForneyLab: SPNonlinearGaussianOutNGP, SPNonlinearGaussianIn1GGP, MNonlinearGaussianGGD

g(x::Float64) = x^2 - 5.0
g(x::Vector{Float64}) = x.^2 .- 5.0
g_inv(x::Float64) = sqrt(x + 5.0)
g_inv(x::Vector{Float64}) = sqrt.(x .+ 5.0)

@testset "SPNonlinearGaussianOutNGP" begin
    @test SPNonlinearGaussianOutNGP <: SumProductRule{NonlinearGaussian}
    @test outboundType(SPNonlinearGaussianOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearGaussianOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, PointMass, m=1.0), g) == Message(Univariate, GaussianMeanVariance, m=2.0000000001164153, v=67.00000000093132)
    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, PointMass, m=1.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.0, v=67.0)
    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.0000000001164153], v=mat(67.00000000093132))
    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(67.0))
end

@testset "SPNonlinearGaussianIn1GGP" begin
    @test SPNonlinearGaussianIn1GGP <: SumProductRule{NonlinearGaussian}
    @test outboundType(SPNonlinearGaussianIn1GGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearGaussianIn1GGP, [Message{Gaussian}, Nothing, Message{PointMass}]) 

    # Without given inverse
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, PointMass, m=1.0), g) == Message(Univariate, GaussianMeanVariance, m=2.499999999868301, v=0.3750000002288054)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, PointMass, m=1.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.5, v=0.375)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.499999999868301], v=mat(0.37500000022152946))
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.5], v=mat(0.375))

    # With given inverse
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, PointMass, m=1.0), g, g_inv) == Message(Univariate, GaussianMeanVariance, m=2.6187538478989154, v=0.14431487180792146)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, PointMass, m=1.0), g, g_inv, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.618033988749895, v=0.14743453366290887)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(1.0)), g, g_inv) == Message(Multivariate, GaussianMeanVariance, m=[2.6187538478989154], v=mat(0.14431487180792146))
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(1.0)), g, g_inv, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.618033988749895], v=mat(0.14743453366290887))
end

@testset "MNonlinearGaussianGGD" begin
    @test MNonlinearGaussianGGD <: MarginalRule{NonlinearGaussian}
    @test isApplicable(MNonlinearGaussianGGD, [Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution])

    @test ruleMNonlinearGaussianGGD(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=1.0), g) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[-0.4444444446018769, 1.388888889228411], w=[1.1111111111424605 -0.22222222229108557; -0.22222222229108557 0.9444444445944998])
    @test ruleMNonlinearGaussianGGD(Message(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0)), g) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[-0.4444444446018769, 1.388888889228411], w=[1.1111111111424605 -0.22222222229108557; -0.22222222229108557 0.9444444445944998])
end

@testset "averageEnergy" begin
    @test averageEnergy(NonlinearGaussian, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)), ProbabilityDistribution(Univariate, PointMass, m=1.0), g) == 10.418938533204672
    @test averageEnergy(NonlinearGaussian, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0)), g) == 10.418938533204672
end

end # module