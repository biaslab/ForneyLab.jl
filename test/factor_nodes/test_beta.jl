module BetaTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeMirroredLogMean, unsafeVar, vague, dims, logPdf
using ForneyLab: SPBetaOutNPP, SPBetaA, SPBetaB, VBBetaOut
using SpecialFunctions: digamma

@testset "Beta ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Univariate, Beta, a=2.0, b=3.0) == ProbabilityDistribution{Univariate, Beta}(Dict(:a=>2.0, :b=>3.0))
    @test_throws Exception ProbabilityDistribution(Multivariate, Beta)
    @test ProbabilityDistribution(Beta, a=2.0, b=3.0) == ProbabilityDistribution{Univariate, Beta}(Dict(:a=>2.0, :b=>3.0))
    @test ProbabilityDistribution(Beta) == ProbabilityDistribution{Univariate, Beta}(Dict(:a=>1.0, :b=>1.0))
    @test Message(Beta) == Message{Beta, Univariate}(ProbabilityDistribution{Univariate, Beta}(Dict(:a=>1.0, :b=>1.0)))
    @test Message(Univariate, Beta) == Message{Beta, Univariate}(ProbabilityDistribution{Univariate, Beta}(Dict(:a=>1.0, :b=>1.0)))
    @test_throws Exception Message(Multivariate, Beta)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Beta, a=2.0, b=2.0)) == 1
end

@testset "vague" begin
    @test vague(Beta) == ProbabilityDistribution(Beta, a=1.0, b=1.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Beta, a=2.0, b=2.0)) == 0.5
    @test unsafeLogMean(ProbabilityDistribution(Beta, a=2.0, b=3.0)) == digamma(2.0) - digamma(5.0)
    @test unsafeMirroredLogMean(ProbabilityDistribution(Beta, a=2.0, b=3.0)) == digamma(3.0) - digamma(5.0)
    @test unsafeVar(ProbabilityDistribution(Beta, a=2.0, b=2.0)) == 0.05
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(Beta, a=2.0, b=2.0),0.3), 0.2311117209633866)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Beta, a=2.0, b=2.0) * ProbabilityDistribution(Beta, a=2.0, b=3.0) == ProbabilityDistribution(Beta, a=3.0, b=4.0)
    @test ProbabilityDistribution(Univariate, Beta, a=1.0, b=2.0) * ProbabilityDistribution(Univariate, PointMass, m=0.2) == ProbabilityDistribution(Univariate, PointMass, m=0.2)
    @test ProbabilityDistribution(Univariate, PointMass, m=0.2) * ProbabilityDistribution(Univariate, Beta, a=1.0, b=2.0) == ProbabilityDistribution(Univariate, PointMass, m=0.2)
    @test_throws Exception ProbabilityDistribution(Univariate, PointMass, m=-1.0) * ProbabilityDistribution(Univariate, Beta, a=1.0, b=2.0)
end


#-------------
# Update rules
#-------------

@testset "SPBetaOutNPP" begin
    @test SPBetaOutNPP <: SumProductRule{Beta}
    @test outboundType(SPBetaOutNPP) == Message{Union{Beta,SampleList}}
    @test isApplicable(SPBetaOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPBetaOutNPP, [Nothing, Message{FactorNode}, Message{PointMass}])
    @test isApplicable(SPBetaOutNPP, [Nothing, Message{PointMass}, Message])

    @test ruleSPBetaOutNPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, Beta, a=2.0, b=3.0)
    d = ruleSPBetaOutNPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, Gamma, a=300.0, b=100.0)).dist
    @test 0.3<mean(d)<0.5
end

@testset "SPBetaA" begin
    @test SPBetaA <: SumProductRule{Beta}
    @test outboundType(SPBetaA) == Message{Function}
    @test isApplicable(SPBetaA, [Message,Nothing,Message])
end

@testset "SPBetaB" begin
    @test SPBetaB <: SumProductRule{Beta}
    @test outboundType(SPBetaB) == Message{Function}
    @test !isApplicable(SPBetaB, [Message,Nothing,Message])
    @test isApplicable(SPBetaB, [Message,Message,Nothing])
end

@testset "VBBetaOut" begin
    @test VBBetaOut <: NaiveVariationalRule{Beta}
    @test outboundType(VBBetaOut) == Message{Beta}
    @test isApplicable(VBBetaOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])

    @test ruleVBBetaOut(nothing, ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Beta, a=2.0, b=3.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test isapprox(differentialEntropy(ProbabilityDistribution(Univariate, Beta, a=2.0, b=3.0)), averageEnergy(Beta, ProbabilityDistribution(Univariate, Beta, a=2.0, b=3.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)))
end

end # module
