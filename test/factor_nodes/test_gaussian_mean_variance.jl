module GaussianMeanVarianceTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPGaussianMeanVarianceOutVPP, SPGaussianMeanVarianceMPVP, SPGaussianMeanVarianceOutVGP, SPGaussianMeanVarianceMGVP, VBGaussianMeanVarianceM, VBGaussianMeanVarianceOut

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
    @test unsafeVar(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeCov(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafePrecision(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.25
    @test unsafeWeightedMean(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.5

    # Multivariate
    @test unsafeMean(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeVar(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [4.0]
    @test unsafeCov(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(4.0)
    @test unsafePrecision(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(0.25)
    @test unsafeWeightedMean(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [0.5]
end


#-------------
# Update rules
#-------------

@testset "SPGaussianMeanVarianceOutVPP" begin
    @test SPGaussianMeanVarianceOutVPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutVPP) == Message{Gaussian}
    @test isApplicable(SPGaussianMeanVarianceOutVPP, [Void, Message{PointMass}, Message{PointMass}]) 
    @test !isApplicable(SPGaussianMeanVarianceOutVPP, [Message{PointMass}, Void, Message{PointMass}]) 

    @test ruleSPGaussianMeanVarianceOutVPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceOutVPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceMPVP" begin
    @test SPGaussianMeanVarianceMPVP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMPVP) == Message{Gaussian}
    @test !isApplicable(SPGaussianMeanVarianceMPVP, [Void, Message{PointMass}, Message{PointMass}]) 
    @test isApplicable(SPGaussianMeanVarianceMPVP, [Message{PointMass}, Void, Message{PointMass}]) 

    @test ruleSPGaussianMeanVarianceMPVP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceMPVP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceOutVGP" begin
    @test SPGaussianMeanVarianceOutVGP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutVGP) == Message{Gaussian}
    @test isApplicable(SPGaussianMeanVarianceOutVGP, [Void, Message{Gaussian}, Message{PointMass}]) 
    @test !isApplicable(SPGaussianMeanVarianceOutVGP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPGaussianMeanVarianceOutVGP(nothing, Message(Univariate, Gaussian, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceOutVGP(nothing, Message(Multivariate, Gaussian, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMeanVarianceMGVP" begin
    @test SPGaussianMeanVarianceMGVP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMGVP) == Message{Gaussian}
    @test !isApplicable(SPGaussianMeanVarianceMGVP, [Void, Message{Gaussian}, Message{PointMass}]) 
    @test isApplicable(SPGaussianMeanVarianceMGVP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPGaussianMeanVarianceMGVP(Message(Univariate, Gaussian, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceMGVP(Message(Multivariate, Gaussian, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))
end

@testset "VBGaussianMeanVarianceM" begin
    @test VBGaussianMeanVarianceM <: NaiveVariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceM) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanVarianceM, [ProbabilityDistribution, Void, ProbabilityDistribution])
    @test !isApplicable(VBGaussianMeanVarianceM, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMeanVarianceM(ProbabilityDistribution(Univariate, Gaussian, m=1.0, v=2.0), nothing, ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceM(ProbabilityDistribution(Multivariate, Gaussian, m=[1.0], v=mat(2.0)), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))
end

@testset "VBGaussianMeanVarianceOut" begin
    @test VBGaussianMeanVarianceOut <: NaiveVariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceOut) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanVarianceOut, [Void, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBGaussianMeanVarianceOut, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMeanVarianceOut(nothing, ProbabilityDistribution(Univariate, Gaussian, m=1.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceOut(nothing, ProbabilityDistribution(Multivariate, Gaussian, m=[1.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=2.0)) == averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=0.0), ProbabilityDistribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=2.0)) == differentialEntropy(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=mat(2.0)))
    @test averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=0.0), ProbabilityDistribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanVariance, ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=mat(2.0)), ProbabilityDistribution(Multivariate, PointMass, m=[0.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)))
end

end #module