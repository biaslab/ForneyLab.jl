module GaussianMeanVarianceTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPGaussianMeanVarianceOutPP, SPGaussianMeanVarianceMPP, SPGaussianMeanVarianceOutGP, SPGaussianMeanVarianceMGP, VBGaussianMeanVarianceM, VBGaussianMeanVarianceOut

#-------------
# Update rules
#-------------

@testset "SPGaussianMeanVarianceOutPP" begin
    @test SPGaussianMeanVarianceOutPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutPP) == Message{Univariate{Gaussian}}
    @test isApplicable(SPGaussianMeanVarianceOutPP, [Void, Message{Univariate{PointMass}}, Message{Univariate{PointMass}}]) 
    @test !isApplicable(SPGaussianMeanVarianceOutPP, [Message{Univariate{PointMass}}, Void, Message{Univariate{PointMass}}]) 

    @test ruleSPGaussianMeanVarianceOutPP(nothing, Message(Univariate(PointMass, m=1.0)), Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=1.0, v=2.0))
end

@testset "SPGaussianMeanVarianceMPP" begin
    @test SPGaussianMeanVarianceMPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMPP) == Message{Univariate{Gaussian}}
    @test !isApplicable(SPGaussianMeanVarianceMPP, [Void, Message{Univariate{PointMass}}, Message{Univariate{PointMass}}]) 
    @test isApplicable(SPGaussianMeanVarianceMPP, [Message{Univariate{PointMass}}, Void, Message{Univariate{PointMass}}]) 

    @test ruleSPGaussianMeanVarianceMPP(Message(Univariate(PointMass, m=1.0)), nothing, Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=1.0, v=2.0))
end

@testset "SPGaussianMeanVarianceOutGP" begin
    @test SPGaussianMeanVarianceOutGP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutGP) == Message{Univariate{Gaussian}}
    @test isApplicable(SPGaussianMeanVarianceOutGP, [Void, Message{Univariate{Gaussian}}, Message{Univariate{PointMass}}]) 
    @test !isApplicable(SPGaussianMeanVarianceOutGP, [Message{Univariate{Gaussian}}, Void, Message{Univariate{PointMass}}]) 

    @test ruleSPGaussianMeanVarianceOutGP(nothing, Message(Univariate(Gaussian, m=1.0, v=1.0)), Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=1.0, v=3.0))
end

@testset "SPGaussianMeanVarianceMGP" begin
    @test SPGaussianMeanVarianceMGP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMGP) == Message{Univariate{Gaussian}}
    @test !isApplicable(SPGaussianMeanVarianceMGP, [Void, Message{Univariate{Gaussian}}, Message{Univariate{PointMass}}]) 
    @test isApplicable(SPGaussianMeanVarianceMGP, [Message{Univariate{Gaussian}}, Void, Message{Univariate{PointMass}}]) 

    @test ruleSPGaussianMeanVarianceMGP(Message(Univariate(Gaussian, m=1.0, v=1.0)), nothing, Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=1.0, v=3.0))
end

@testset "VBGaussianMeanVarianceM" begin
    @test VBGaussianMeanVarianceM <: VariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceM) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMeanVarianceM, [Univariate, Void, Univariate])
    @test !isApplicable(VBGaussianMeanVarianceM, [Univariate, Univariate, Void]) 

    @test ruleVBGaussianMeanVarianceM(Univariate(Gaussian, m=1.0, v=2.0), nothing, Univariate(PointMass, m=3.0)) == Message(Univariate(Gaussian, m=1.0, v=3.0))
end

@testset "VBGaussianMeanVarianceOut" begin
    @test VBGaussianMeanVarianceOut <: VariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceOut) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMeanVarianceOut, [Void, Univariate, Univariate])
    @test !isApplicable(VBGaussianMeanVarianceOut, [Univariate, Univariate, Void]) 

    @test ruleVBGaussianMeanVarianceOut(nothing, Univariate(Gaussian, m=1.0, v=2.0), Univariate(PointMass, m=3.0)) == Message(Univariate(Gaussian, m=1.0, v=3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Gaussian, m=0.0, v=2.0)) == averageEnergy(GaussianMeanVariance, Univariate(Gaussian, m=0.0, v=2.0), Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0))
    @test differentialEntropy(Univariate(Gaussian, m=0.0, v=2.0)) == differentialEntropy(Multivariate(Gaussian, m=[0.0], v=[2.0].'))
    @test averageEnergy(GaussianMeanVariance, Univariate(Gaussian, m=0.0, v=2.0), Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0)) == averageEnergy(GaussianMeanVariance, Multivariate(Gaussian, m=[0.0], v=[2.0].'), Multivariate(PointMass, m=[0.0]), MatrixVariate(PointMass, m=[2.0].'))
end

end #module