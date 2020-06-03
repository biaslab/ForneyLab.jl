module GaussianMeanVarianceTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianMeanVarianceOutNPP, SPGaussianMeanVarianceMPNP, SPGaussianMeanVarianceOutNGP, SPGaussianMeanVarianceMGNP, SPGaussianMeanVarianceVGGN, SPGaussianMeanVarianceVPGN, VBGaussianMeanVarianceM, VBGaussianMeanVarianceOut
using LinearAlgebra: det, diag

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=ones(2), v=diageye(2))) == 2
end

@testset "vague" begin
    @test vague(GaussianMeanVariance) == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
    @test vague(GaussianMeanVariance, 2) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=huge*eye(2))
    @test vague(GaussianMeanVariance, (2,)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=huge*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0))
    @test !isProper(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=-1.0))

    # Multivariate
    @test isProper(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)))
    @test isProper(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=ones(2), v=diageye(2)))
    @test !isProper(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    @test ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)

    # Multivariate
    @test ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    @test ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 2.0
    @test unsafeMode(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 2.0
    @test unsafeVar(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeCov(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeMeanCov(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (2.0, 4.0)
    @test unsafePrecision(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.25
    @test unsafeWeightedMean(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.5
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (0.5, 0.25)

    # Multivariate
    @test unsafeMean(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeMode(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeVar(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [4.0]
    @test unsafeCov(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(4.0)
    @test unsafeMeanCov(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([2.0], mat(4.0))
    @test unsafePrecision(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(0.25)
    @test unsafeWeightedMean(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [0.5]
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([0.5], mat(0.25))
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=0.5), 1.0), -0.5723649429247)
    @test isapprox(logPdf(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0, 1.0], v=[0.5 0.0; 0.0 0.5]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=8.0, w=4.0)) == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=0.25)
    @test convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=0.25)
    @test convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[8.0], w=mat(4.0))) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.25))
    @test convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.25))
end


#-------------
# Update rules
#-------------

@testset "SPGaussianMeanVarianceOutNPP" begin
    @test SPGaussianMeanVarianceOutNPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutNPP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanVarianceOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutNPP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceMPNP" begin
    @test SPGaussianMeanVarianceMPNP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMPNP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanVarianceMPNP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanVarianceMPNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceMPNP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceMPNP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceOutNGP" begin
    @test SPGaussianMeanVarianceOutNGP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanVarianceOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutNGP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMeanVarianceMGNP" begin
    @test SPGaussianMeanVarianceMGNP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMGNP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanVarianceMGNP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanVarianceMGNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceMGNP(Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceMGNP(Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMeanVarianceVGGN" begin
    @test SPGaussianMeanVarianceVGGN <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVGGN) == Message{Function}
    @test !isApplicable(SPGaussianMeanVarianceVGGN, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test isApplicable(SPGaussianMeanVarianceVGGN, [Message{Gaussian}, Message{Gaussian}, Nothing])

    msg = ruleSPGaussianMeanVarianceVGGN(Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing)
    @test isa(msg, Message{Function, Univariate})
    @test msg.dist.params[:log_pdf](1.0) == -0.5*log(2.0 + 4.0 + 1.0) - 1/(2*1.0)*(1.0 - 3.0)^2
end

@testset "SPGaussianMeanVarianceVPGN" begin
    @test SPGaussianMeanVarianceVPGN <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVPGN) == Message{Function}
    @test isApplicable(SPGaussianMeanVarianceVPGN, [Message{PointMass}, Message{Gaussian}, Nothing])

    msg = ruleSPGaussianMeanVarianceVPGN(Message(Univariate, PointMass, m=1.0), Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing)
    @test isa(msg, Message{Function, Univariate})
    @test msg.dist.params[:log_pdf](1.0) == -0.5*log(4.0 + 1.0) - 1/(2*1.0)*(1.0 - 3.0)^2
end

@testset "VBGaussianMeanVarianceM" begin
    @test VBGaussianMeanVarianceM <: NaiveVariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceM) == Message{GaussianMeanVariance}
    @test isApplicable(VBGaussianMeanVarianceM, [ProbabilityDistribution, Nothing, ProbabilityDistribution])
    @test !isApplicable(VBGaussianMeanVarianceM, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBGaussianMeanVarianceM(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing, ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceM(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "VBGaussianMeanVarianceOut" begin
    @test VBGaussianMeanVarianceOut <: NaiveVariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceOut) == Message{GaussianMeanVariance}
    @test isApplicable(VBGaussianMeanVarianceOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBGaussianMeanVarianceOut, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBGaussianMeanVarianceOut(nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceOut(nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0)) == averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=0.0), ProbabilityDistribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0)) == differentialEntropy(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(2.0)))
    @test averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=0.0), ProbabilityDistribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(2.0)), ProbabilityDistribution(Multivariate, PointMass, m=[0.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)))
    @test averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), ProbabilityDistribution(Univariate, PointMass, m=0.5)) == averageEnergy(GaussianMeanPrecision, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), ProbabilityDistribution(Univariate, PointMass, m=2.0))
end

end #module
