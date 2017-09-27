module GammaTest

using Base.Test
using ForneyLab
import ForneyLab: prod!, unsafeMean, unsafeVar, outboundType, isApplicable, VBGamma3

@testset "prod!" begin
    @test ProbabilityDistribution(Gamma, a=1.0, b=2.0) * ProbabilityDistribution(Gamma, a=3.0, b=4.0) == ProbabilityDistribution(Gamma, a=3.0, b=6.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Gamma, a=1.0, b=2.0)) == 0.5
    @test unsafeVar(ProbabilityDistribution(Gamma, a=1.0, b=2.0)) == 0.25
end

#-------------
# Update rules
#-------------

@testset "VBGamma3" begin
    @test VBGamma3 <: VariationalRule{Gamma}
    @test outboundType(VBGamma3) == Message{Gamma}
    @test isApplicable(VBGamma3, 3) 
    @test !isApplicable(VBGamma3, 2) 

    @test ruleVBGamma3(ProbabilityDistribution(PointMass, m=1.5), ProbabilityDistribution(PointMass, m=3.0), nothing) == Message(Gamma, a=1.5, b=3.0)
end


end #module