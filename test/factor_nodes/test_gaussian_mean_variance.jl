module GaussianMeanVarianceTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPGaussianMeanVariancePPV, SPGaussianMeanVarianceVPP, SPGaussianMeanVarianceGPV, SPGaussianMeanVarianceVPG, VBGaussianMeanVariance3


#-------------
# Update rules
#-------------

@testset "SPGaussianMeanVariancePPV" begin
    @test SPGaussianMeanVariancePPV <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVariancePPV) == Message{Gaussian}
    @test isApplicable(SPGaussianMeanVariancePPV, [Message{PointMass}, Message{PointMass}, Void]) 
    @test !isApplicable(SPGaussianMeanVariancePPV, [Void, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPGaussianMeanVariancePPV(Message(PointMass, m=1.0), Message(PointMass, m=2.0), nothing) == Message(Gaussian, m=1.0, v=2.0)
end

@testset "SPGaussianMeanVarianceVPP" begin
    @test SPGaussianMeanVarianceVPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVPP) == Message{Gaussian}
    @test !isApplicable(SPGaussianMeanVarianceVPP, [Message{PointMass}, Message{PointMass}, Void]) 
    @test isApplicable(SPGaussianMeanVarianceVPP, [Void, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPGaussianMeanVarianceVPP(nothing, Message(PointMass, m=2.0), Message(PointMass, m=1.0)) == Message(Gaussian, m=1.0, v=2.0)
end

@testset "SPGaussianMeanVarianceGPV" begin
    @test SPGaussianMeanVarianceGPV <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceGPV) == Message{Gaussian}
    @test isApplicable(SPGaussianMeanVarianceGPV, [Message{Gaussian}, Message{PointMass}, Void]) 
    @test !isApplicable(SPGaussianMeanVarianceGPV, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPGaussianMeanVarianceGPV(Message(Gaussian, m=1.0, v=1.0), Message(PointMass, m=2.0), nothing) == Message(Gaussian, m=1.0, v=3.0)
end

@testset "SPGaussianMeanVarianceVPG" begin
    @test SPGaussianMeanVarianceVPG <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVPG) == Message{Gaussian}
    @test !isApplicable(SPGaussianMeanVarianceVPG, [Message{Gaussian}, Message{PointMass}, Void]) 
    @test isApplicable(SPGaussianMeanVarianceVPG, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPGaussianMeanVarianceVPG(nothing, Message(PointMass, m=2.0), Message(Gaussian, m=1.0, v=1.0)) == Message(Gaussian, m=1.0, v=3.0)
end

@testset "VBGaussianMeanVariance3" begin
    @test VBGaussianMeanVariance3 <: VariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVariance3) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanVariance3, GaussianMeanVariance, 3) 
    @test !isApplicable(VBGaussianMeanVariance3, GaussianMeanVariance, 2) 

    @test ruleVBGaussianMeanVariance3(ProbabilityDistribution(Gaussian, m=1.0, v=2.0), ProbabilityDistribution(PointMass, m=3.0), nothing) == Message(Gaussian, m=1.0, v=3.0)
end

end #module