module GaussianMeanVarianceTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
import ForneyLab: SPGaussianMeanVarianceOutVPP, SPGaussianMeanVarianceMPVP, SPGaussianMeanVarianceOutVGP, SPGaussianMeanVarianceMGVP, VBGaussianMeanVarianceM, VBGaussianMeanVarianceOut
import PDMats: PDMat, PDiagMat
import LinearAlgebra: I, det, diag, logdet

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=ones(2), v=diageye(2))) == 2
end

@testset "vague" begin
    @test vague(GaussianMeanVariance) == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
    @test vague(GaussianMeanVariance, 2) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=huge*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0))
    @test !isProper(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=-1.0))

    # Multivariate
    @test isProper(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)))
    @test isProper(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=ones(2), v=diageye(2)))
    # @test !isProper(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(-1.0)))
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
    @test unsafeVar(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeCov(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeMeanCov(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (2.0, 4.0)
    @test unsafePrecision(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.25
    @test unsafeWeightedMean(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.5
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (0.5, 0.25)

    # Multivariate
    @test unsafeMean(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeVar(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [4.0]
    @test unsafeCov(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == PDMat(mat(4.0))
    @test unsafeMeanCov(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([2.0], PDMat(mat(4.0)))
    @test unsafePrecision(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(0.25)
    @test unsafeWeightedMean(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [0.5]
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([0.5], mat(0.25))
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

@testset "SPGaussianMeanVarianceOutVPP" begin
    @test SPGaussianMeanVarianceOutVPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutVPP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanVarianceOutVPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutVPP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutVPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceOutVPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceMPVP" begin
    @test SPGaussianMeanVarianceMPVP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMPVP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanVarianceMPVP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanVarianceMPVP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceMPVP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceMPVP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceOutVGP" begin
    @test SPGaussianMeanVarianceOutVGP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutVGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanVarianceOutVGP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutVGP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutVGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceOutVGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMeanVarianceMGVP" begin
    @test SPGaussianMeanVarianceMGVP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMGVP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanVarianceMGVP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanVarianceMGVP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceMGVP(Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceMGVP(Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
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
end

end #module