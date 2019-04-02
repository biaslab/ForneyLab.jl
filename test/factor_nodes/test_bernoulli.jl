module BernoulliTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, vague, dims
import ForneyLab: SPBernoulliOutVP, SPBernoulliIn1PV, SPBernoulliOutVB, VBBernoulliOut, VBBernoulliIn1

@testset "Bernoulli ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Univariate, Bernoulli, p=0.2) == ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>0.2))
    @test_throws Exception ProbabilityDistribution(Multivariate, Bernoulli)
    @test ProbabilityDistribution(Bernoulli, p=0.2) == ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>0.2))
    @test ProbabilityDistribution(Bernoulli) == ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>0.5))
    @test Message(Bernoulli) == Message{Bernoulli, Univariate}(ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>0.5)))
    @test Message(Univariate, Bernoulli) == Message{Bernoulli, Univariate}(ProbabilityDistribution{Univariate, Bernoulli}(Dict(:p=>0.5)))
    @test_throws Exception Message(Multivariate, Bernoulli)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Bernoulli, p=0.5)) == 1
end

@testset "vague" begin
    @test vague(Bernoulli) == ProbabilityDistribution(Bernoulli, p=0.5)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Bernoulli, p=0.2)) == 0.2
    @test unsafeVar(ProbabilityDistribution(Bernoulli, p=0.5)) == 0.25
end

@testset "prod!" begin
    @test ProbabilityDistribution(Bernoulli, p=0.2) * ProbabilityDistribution(Bernoulli, p=0.8) == ProbabilityDistribution(Bernoulli, p=0.5000000000000001)
    @test_throws Exception ProbabilityDistribution(Bernoulli, p=0.0) * ProbabilityDistribution(Bernoulli, p=1.0)
end

#-------------
# Update rules
#-------------

@testset "SPBernoulliOutVP" begin
    @test SPBernoulliOutVP <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliOutVP) == Message{Bernoulli}
    @test isApplicable(SPBernoulliOutVP, [Nothing, Message{PointMass}])

    @test ruleSPBernoulliOutVP(nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, Bernoulli, p=1.0)
end

@testset "SPBernoulliIn1PV" begin
    @test SPBernoulliIn1PV <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliIn1PV) == Message{Beta}
    @test isApplicable(SPBernoulliIn1PV, [Message{PointMass}, Nothing])

    @test ruleSPBernoulliIn1PV(Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, Beta, a=2.0, b=1.0)
end

@testset "SPBernoulliOutVB" begin
    @test SPBernoulliOutVB <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliOutVB) == Message{Bernoulli}
    @test isApplicable(SPBernoulliOutVB, [Nothing, Message{Beta}])

    @test ruleSPBernoulliOutVB(nothing, Message(Univariate, Beta, a=1.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.5)
end

@testset "VBBernoulliOut" begin
    @test VBBernoulliOut <: NaiveVariationalRule{Bernoulli}
    @test outboundType(VBBernoulliOut) == Message{Bernoulli}
    @test isApplicable(VBBernoulliOut, [Nothing, ProbabilityDistribution])
    @test !isApplicable(VBBernoulliOut, [ProbabilityDistribution, Nothing])

    @test ruleVBBernoulliOut(nothing, ProbabilityDistribution(Univariate, PointMass, m=0.2)) == Message(Univariate, Bernoulli, p=0.2)
    @test ruleVBBernoulliOut(nothing, ProbabilityDistribution(Univariate, Beta, a=1.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.5)
end

@testset "VBBernoulliIn1" begin
    @test VBBernoulliIn1 <: NaiveVariationalRule{Bernoulli}
    @test outboundType(VBBernoulliIn1) == Message{Beta}
    @test isApplicable(VBBernoulliIn1, [ProbabilityDistribution, Nothing])

    @test ruleVBBernoulliIn1(ProbabilityDistribution(Univariate, Bernoulli, p=0.2), nothing) == Message(Univariate, Beta, a=1.2, b=1.8)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, Bernoulli, p=0.25)) == averageEnergy(Bernoulli, ProbabilityDistribution(Univariate, Bernoulli, p=0.25), ProbabilityDistribution(Univariate, PointMass, m=0.25))
end

end # module
