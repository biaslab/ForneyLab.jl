module GammaTest

using Test
using ForneyLab
import ForneyLab: prod!, unsafeMean, unsafeVar, outboundType, isApplicable, dims
import ForneyLab: SPGammaOutVPP, VBGammaOut

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=1.0)) == 1
end

@testset "vague" begin
    @test vague(Gamma) == ProbabilityDistribution(Univariate, Gamma, a=1.0, b=tiny)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0) * ProbabilityDistribution(Univariate, Gamma, a=3.0, b=4.0) == ProbabilityDistribution(Univariate, Gamma, a=3.0, b=6.0)
    @test ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0) * ProbabilityDistribution(Univariate, PointMass, m=1.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)
    @test ProbabilityDistribution(Univariate, PointMass, m=1.0) * ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)
    @test_throws Exception ProbabilityDistribution(Univariate, PointMass, m=-1.0) * ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == 0.5
    @test unsafeVar(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == 0.25
end


#-------------
# Update rules
#-------------

@testset "SPGammaOutVPP" begin
    @test SPGammaOutVPP <: SumProductRule{Gamma}
    @test outboundType(SPGammaOutVPP) == Message{Gamma}
    @test isApplicable(SPGammaOutVPP, [Nothing, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPGammaOutVPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gamma, a=1.0, b=2.0)
end

@testset "VBGammaOut" begin
    @test VBGammaOut <: NaiveVariationalRule{Gamma}
    @test outboundType(VBGammaOut) == Message{Gamma}
    @test isApplicable(VBGammaOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGammaOut, [ProbabilityDistribution, ProbabilityDistribution, Nothing]) 

    @test ruleVBGammaOut(nothing, ProbabilityDistribution(Univariate, PointMass, m=1.5), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gamma, a=1.5, b=3.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == averageEnergy(Gamma, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, PointMass, m=1.0), ProbabilityDistribution(Univariate, PointMass, m=2.0))
end

end #module