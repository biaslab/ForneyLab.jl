module GammaTest

using Base.Test
using ForneyLab
import ForneyLab: prod!, unsafeMean, unsafeVar, outboundType, isApplicable
import ForneyLab: VBGammaOut

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
    @test outboundType(VBGammaOut) == Message{Univariate{Gamma}}
    @test isApplicable(VBGammaOut, [Void, Univariate, Univariate]) 
    @test !isApplicable(VBGammaOut, [Univariate, Univariate, Void]) 

    @test ruleVBGammaOut(nothing, Univariate(PointMass, m=1.5), Univariate(PointMass, m=3.0)) == Message(Univariate(Gamma, a=1.5, b=3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Gamma, a=1.0, b=2.0)) == averageEnergy(Gamma, Univariate(Gamma, a=1.0, b=2.0), Univariate(PointMass, m=1.0), Univariate(PointMass, m=2.0))
end

end #module