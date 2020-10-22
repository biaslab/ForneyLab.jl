module BetaTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeMirroredLogMean, unsafeVar, vague, dims, logPdf
using ForneyLab: SPBetaOutNMM, SPBetaMNM, SPBetaMMN, SPBetaOutMCNMM, VBBetaOut, VBBetaA, VBBetaB
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

@testset "SPBetaOutNMM" begin
    @test SPBetaOutNMM <: SumProductRule{Beta}
    @test outboundType(SPBetaOutNMM) == Message{Beta}
    @test isApplicable(SPBetaOutNMM, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPBetaOutNMM, [Nothing, Message{FactorNode}, Message{PointMass}])
    @test !isApplicable(SPBetaOutNMM, [Nothing, Message{PointMass}, Message{FactorNode}])

    @test ruleSPBetaOutNMM(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, Beta, a=2.0, b=3.0)
end

@testset "SPBetaOutMCNMM" begin
    @test SPBetaOutMCNMM <: SumProductRule{Beta}
    @test outboundType(SPBetaOutMCNMM) == Message{SampleList}
    @test !isApplicable(SPBetaOutMCNMM, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPBetaOutMCNMM, [Nothing, Message{FactorNode}, Message{PointMass}])
    @test isApplicable(SPBetaOutMCNMM, [Nothing, Message{PointMass}, Message{FactorNode}])

    d = ruleSPBetaOutMCNMM(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, Gamma, a=300.0, b=100.0)).dist
    @test 0.3<mean(d)<0.5
end

@testset "SPBetaMNM" begin
    @test SPBetaMNM <: SumProductRule{Beta}
    @test outboundType(SPBetaMNM) == Message{Function}
    @test isApplicable(SPBetaMNM, [Message,Nothing,Message])
end

@testset "SPBetaMMN" begin
    @test SPBetaMMN <: SumProductRule{Beta}
    @test outboundType(SPBetaMMN) == Message{Function}
    @test !isApplicable(SPBetaMMN, [Message,Nothing,Message])
    @test isApplicable(SPBetaMMN, [Message,Message,Nothing])
end

@testset "VBBetaOut" begin
    @test VBBetaOut <: NaiveVariationalRule{Beta}
    @test outboundType(VBBetaOut) == Message{Beta}
    @test isApplicable(VBBetaOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])

    @test ruleVBBetaOut(nothing, ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Beta, a=2.0, b=3.0)
end

@testset "VBBetaA" begin
    @test VBBetaA <: NaiveVariationalRule{Beta}
    @test outboundType(VBBetaA) == Message{Function}
    @test isApplicable(VBBetaA, [ProbabilityDistribution, Nothing, ProbabilityDistribution])
end

@testset "VBBetaB" begin
    @test VBBetaB <: NaiveVariationalRule{Beta}
    @test outboundType(VBBetaB) == Message{Function}
    @test isApplicable(VBBetaB, [ProbabilityDistribution, ProbabilityDistribution, Nothing])
end

@testset "averageEnergy and differentialEntropy" begin
    @test isapprox(differentialEntropy(ProbabilityDistribution(Univariate, Beta, a=2.0, b=3.0)), averageEnergy(Beta, ProbabilityDistribution(Univariate, Beta, a=2.0, b=3.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)))
end

end # module
