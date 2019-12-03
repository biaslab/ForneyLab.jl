module BetaTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeMirroredLogMean, unsafeVar, vague, dims, logPdf
import ForneyLab: SPBetaOutNPP, VBBetaOut
import SpecialFunctions: digamma

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
    @test outboundType(SPBetaOutNPP) == Message{Beta}
    @test isApplicable(SPBetaOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPBetaOutNPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, Beta, a=2.0, b=3.0)
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
