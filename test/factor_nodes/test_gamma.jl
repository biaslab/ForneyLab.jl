module GammaTest

using Base.Test
using ForneyLab
import ForneyLab: prod!, unsafeMean, unsafeVar, outboundType, isApplicable, dims
import ForneyLab: VBGammaOut

@testset "dims" begin
    @test dims(Univariate(Gamma, a=1.0, b=1.0)) == 1
end

@testset "vague" begin
    @test vague(Gamma) == Univariate(Gamma, a=1.0, b=tiny)
end

@testset "prod!" begin
    @test Univariate(Gamma, a=1.0, b=2.0) * Univariate(Gamma, a=3.0, b=4.0) == Univariate(Gamma, a=3.0, b=6.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Univariate(Gamma, a=1.0, b=2.0)) == 0.5
    @test unsafeVar(Univariate(Gamma, a=1.0, b=2.0)) == 0.25
end

#-------------
# Update rules
#-------------

@testset "VBGammaOut" begin
    @test VBGammaOut <: VariationalRule{Gamma}
    @test outboundType(VBGammaOut) == Message{AbstractGamma}
    @test isApplicable(VBGammaOut, [Void, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGammaOut, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGammaOut(nothing, Univariate(PointMass, m=1.5), Univariate(PointMass, m=3.0)) == Message(Univariate(Gamma, a=1.5, b=3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Gamma, a=1.0, b=2.0)) == averageEnergy(Gamma, Univariate(Gamma, a=1.0, b=2.0), Univariate(PointMass, m=1.0), Univariate(PointMass, m=2.0))
end

end #module