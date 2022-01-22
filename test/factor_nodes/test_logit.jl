module LogitTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: VBLogitOut, VBLogitIn1, VBLogitXi


#-------------
# Update rules
#-------------

@testset "VBLogitOut" begin
    @test VBLogitOut <: NaiveVariationalRule{Logit}
    @test outboundType(VBLogitOut) == Message{Bernoulli}
    @test isApplicable(VBLogitOut, [Nothing, Distribution, Distribution]) 

    @test ruleVBLogitOut(nothing, Distribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Bernoulli, p=1/(1+exp(-2.0)))
end

@testset "VBLogitIn1" begin
    @test VBLogitIn1 <: NaiveVariationalRule{Logit}
    @test outboundType(VBLogitIn1) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(VBLogitIn1, [Distribution, Nothing, Distribution]) 

    @test ruleVBLogitIn1(Distribution(Univariate, Bernoulli, p=0.8), nothing, Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.3, w=(1/(1+exp(-3.0)) - 0.5)/3.0)
end

@testset "VBLogitXi" begin
    @test VBLogitXi <: NaiveVariationalRule{Logit}
    @test outboundType(VBLogitXi) == Message{Function}
    @test isApplicable(VBLogitXi, [Distribution, Distribution, Nothing]) 

    @test ruleVBLogitXi(Distribution(Univariate, Bernoulli, p=0.8), Distribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), nothing) == Message(Univariate, Function, mode=sqrt(5.0))
end

@testset "averageEnergy" begin
    @test averageEnergy(Logit, Distribution(Univariate, Bernoulli, p=0.8), Distribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Distribution(Univariate, PointMass, m=3.0)) == (1/(1+exp(-3.0)) - 0.5)/6.0*(5.0 - 9.0) + 2.5 + log(1+exp(-3.0)) - 1.6 
end

end # module