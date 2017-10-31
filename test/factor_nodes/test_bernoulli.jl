module BernoulliTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, vague, dims
import ForneyLab: SPBernoulliOutVP, VBBernoulliOut

@testset "dims" begin
    @test dims(Univariate(Bernoulli, p=0.5)) == 1
end

@testset "vague" begin
    @test vague(Bernoulli) == Univariate(Bernoulli, p=0.5)
end

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

@testset "SPBernoulliOutVP" begin
    @test SPBernoulliOutVP <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliOutVP) == Message{Bernoulli}
    @test isApplicable(SPBernoulliOutVP, [Void, Message{PointMass}]) 

    @test ruleSPBernoulliOutVP(nothing, Message(Univariate(PointMass, m=0.2))) == Message(Univariate(Bernoulli, p=0.2))
end

@testset "VBBernoulliOut" begin
    @test VBBernoulliOut <: VariationalRule{Bernoulli}
    @test outboundType(VBBernoulliOut) == Message{Bernoulli}
    @test isApplicable(VBBernoulliOut, [Void, ProbabilityDistribution])
    @test !isApplicable(VBBernoulliOut, [ProbabilityDistribution, Void])

    @test ruleVBBernoulliOut(nothing, Univariate(PointMass, m=0.2)) == Message(Univariate(Bernoulli, p=0.2))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Bernoulli, p=0.25)) == averageEnergy(Bernoulli, Univariate(Bernoulli, p=0.25), Univariate(PointMass, m=0.25))
end

end # module