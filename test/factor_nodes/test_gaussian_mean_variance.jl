module GaussianMeanVarianceTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPGaussianMeanVarianceOutVPP, SPGaussianMeanVarianceMPVP, SPGaussianMeanVarianceOutVGP, SPGaussianMeanVarianceMGVP, VBGaussianMeanVarianceM, VBGaussianMeanVarianceOut

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
    @test VBGaussianMeanVarianceM <: VariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceM) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanVarianceM, [ProbabilityDistribution, Void, ProbabilityDistribution])
    @test !isApplicable(VBGaussianMeanVarianceM, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMeanVarianceM(ProbabilityDistribution(Univariate, Gaussian, m=1.0, v=2.0), nothing, ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceM(ProbabilityDistribution(Multivariate, Gaussian, m=[1.0], v=mat(2.0)), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))
end

@testset "VBGaussianMeanVarianceOut" begin
    @test VBGaussianMeanVarianceOut <: VariationalRule{GaussianMeanVariance}
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