module BernoulliTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar
import ForneyLab: SPBernoulliOutP, VBBernoulliOut

@testset "unsafe mean and variance" begin
    @test unsafeMean(Univariate(Bernoulli, p=0.2)) == 0.2
    @test unsafeVar(Univariate(Bernoulli, p=0.5)) == 0.25
end

@testset "prod!" begin
    @test Univariate(Bernoulli, p=0.2) * Univariate(Bernoulli, p=0.8) == Univariate(Bernoulli, p=0.5000000000000001)
    @test_throws Exception Univariate(Bernoulli, p=0.0) * Univariate(Bernoulli, p=1.0)
end

#-------------
# Update rules
#-------------

@testset "SPBernoulliOutP" begin
    @test SPBernoulliOutP <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliOutP) == Message{Univariate{Bernoulli}}
    @test isApplicable(SPBernoulliOutP, [Void, Message{Univariate{PointMass}}]) 

    @test ruleSPBernoulliOutP(nothing, Message(Univariate(PointMass, m=0.2))) == Message(Univariate(Bernoulli, p=0.2))
end

@testset "VBBernoulliOut" begin
    @test VBBernoulliOut <: VariationalRule{Bernoulli}
    @test outboundType(VBBernoulliOut) == Message{Univariate{Bernoulli}}
    @test isApplicable(VBBernoulliOut, [Void, Univariate])
    @test !isApplicable(VBBernoulliOut, [Univariate, Void])

    @test ruleVBBernoulliOut(nothing, Univariate(PointMass, m=0.2)) == Message(Univariate(Bernoulli, p=0.2))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Bernoulli, p=0.25)) == averageEnergy(Bernoulli, Univariate(Bernoulli, p=0.25), Univariate(PointMass, m=0.25))
end

end # module