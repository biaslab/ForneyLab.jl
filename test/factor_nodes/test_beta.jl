module BetaTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeMirroredLogMean, unsafeVar, vague, dims, logPdf, naturalParams, standardDistribution
using ForneyLab: SPBetaOutNPP, SPBetaAMNM, SPBetaBMMN, SPBetaOutNMM, VBBetaOut, VBBetaA, VBBetaB
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
    @test dims(ProbabilityDistribution(Beta, a=2.0, b=2.0)) == ()
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

@testset "natural parameters" begin
    d = ProbabilityDistribution(Univariate, Beta, a=2.0, b=3.0)
    η = naturalParams(d)
    s = standardDistribution(Univariate, Beta, η=η)
    @test d.params[:a] == s.params[:a] # Test conversion consistency
    @test d.params[:b] == s.params[:b]

    x = [0.1, 0.6, 0.9]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, Beta, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
end


#-------------
# Update rules
#-------------

@testset "SPBetaOutNPP" begin
    @test SPBetaOutNPP <: SumProductRule{Beta}
    @test outboundType(SPBetaOutNPP) == Message{Beta}
    @test isApplicable(SPBetaOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPBetaOutNPP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPBetaOutNPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, Beta, a=2.0, b=3.0)
end

@testset "SPBetaOutNMM" begin
    @test SPBetaOutNMM <: SumProductRule{Beta}
    @test outboundType(SPBetaOutNMM) == Message{SampleList}
    @test !isApplicable(SPBetaOutNMM, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPBetaOutNMM, [Nothing, Message{FactorNode}, Message{PointMass}])
    @test isApplicable(SPBetaOutNMM, [Nothing, Message{PointMass}, Message{FactorNode}])

    msg = ruleSPBetaOutNMM(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, Gamma, a=300.0, b=100.0))
    @test length(msg.dist.params[:s]) == length(msg.dist.params[:w])
end

@testset "SPBetaAMNM" begin
    @test SPBetaAMNM <: SumProductRule{Beta}
    @test outboundType(SPBetaAMNM) == Message{Function}
    @test isApplicable(SPBetaAMNM, [Message, Nothing, Message])
end

@testset "SPBetaBMMN" begin
    @test SPBetaBMMN <: SumProductRule{Beta}
    @test outboundType(SPBetaBMMN) == Message{Function}
    @test !isApplicable(SPBetaBMMN, [Message, Nothing, Message])
    @test isApplicable(SPBetaBMMN, [Message, Message, Nothing])
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
