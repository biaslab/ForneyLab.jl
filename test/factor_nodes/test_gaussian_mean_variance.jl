module GaussianMeanVarianceTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPGaussianMeanVariancePPV, SPGaussianMeanVarianceVPP


#-------------
# Update rules
#-------------

@testset "SPGaussianMeanVariancePPV" begin
    @test SPGaussianMeanVariancePPV <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVariancePPV) == Message{Gaussian}
    @test isApplicable(SPGaussianMeanVariancePPV, [Message{PointMass}, Message{PointMass}, Void]) 
    @test !isApplicable(SPGaussianMeanVariancePPV, [Void, Message{PointMass}, Message{PointMass}]) 
end

@testset "SPGaussianMeanVarianceVPP" begin
    @test SPGaussianMeanVarianceVPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVPP) == Message{Gaussian}
    @test !isApplicable(SPGaussianMeanVarianceVPP, [Message{PointMass}, Message{PointMass}, Void]) 
    @test isApplicable(SPGaussianMeanVarianceVPP, [Void, Message{PointMass}, Message{PointMass}]) 
end

# TODO: Add more tests
@testset "ruleSPGaussianMeanVariancePPV" begin
    @test ruleSPGaussianMeanVariancePPV(Message(PointMass, m=1.0), Message(PointMass, m=2.0), nothing) == Message(Gaussian, m=1.0, v=2.0)
end

@testset "ruleSPGaussianMeanVarianceVPP" begin
    @test ruleSPGaussianMeanVarianceVPP(nothing, Message(PointMass, m=2.0), Message(PointMass, m=1.0)) == Message(Gaussian, m=1.0, v=2.0)
end

end #module