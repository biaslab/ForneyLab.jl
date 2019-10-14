module LogNormalTest

using Test
using ForneyLab
import ForneyLab: prod!, unsafeMean, unsafeLogMean, unsafeVar, unsafeLogVar, unsafeCov, unsafeLogCov, outboundType, isApplicable, dims
import ForneyLab: SPLogNormalOutNPP, VBLogNormalOut

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=1.0)) == 1
end

@testset "vague" begin
    @test vague(LogNormal) == ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=huge)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == exp(1.0 + 0.5*2.0)
    @test unsafeLogMean(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == 1.0
    @test unsafeVar(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == exp(2.0*1.0 + 2.0)*(exp(2.0) - 1.0)
    @test unsafeLogVar(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == 2.0
    @test unsafeCov(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == exp(2.0*1.0 + 2.0)*(exp(2.0) - 1.0)
    @test unsafeLogCov(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == 2.0
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(Univariate, LogNormal, m=1.2, s=0.5),2), -1.522411904058978)
end

@testset "Gamma approximatons to LogNormal" begin
    @test ForneyLab.laplace(Gamma, ProbabilityDistribution(Univariate, LogNormal, m=0.0, s=2.0)) == ProbabilityDistribution(Univariate, Gamma, a=0.5, b=0.5)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0) * ProbabilityDistribution(Univariate, PointMass, m=1.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)
    @test_throws Exception ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0) * ProbabilityDistribution(Univariate, PointMass, m=-1.0)
    @test ProbabilityDistribution(Univariate, LogNormal, m=0.0, s=2.0) * ProbabilityDistribution(Univariate, Gamma, a=3.0, b=4.0) == ProbabilityDistribution(Univariate, Gamma, a=2.5, b=4.5)
end


#-------------
# Update rules
#-------------

@testset "SPLogNormalOutNPP" begin
    @test SPLogNormalOutNPP <: SumProductRule{LogNormal}
    @test outboundType(SPLogNormalOutNPP) == Message{LogNormal}
    @test isApplicable(SPLogNormalOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPLogNormalOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, LogNormal, m=1.0, s=2.0)
end

@testset "VBLogNormalOut" begin
    @test VBLogNormalOut <: NaiveVariationalRule{LogNormal}
    @test outboundType(VBLogNormalOut) == Message{LogNormal}
    @test isApplicable(VBLogNormalOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBLogNormalOut, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBLogNormalOut(nothing, ProbabilityDistribution(Univariate, PointMass, m=1.5), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, LogNormal, m=1.5, s=3.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0)) == averageEnergy(LogNormal, ProbabilityDistribution(Univariate, LogNormal, m=1.0, s=2.0), ProbabilityDistribution(Univariate, PointMass, m=1.0), ProbabilityDistribution(Univariate, PointMass, m=2.0))
end

end #module
