module LogNormalTest

using Test
using ForneyLab
using ForneyLab: prod!, unsafeMean, unsafeLogMean, unsafeVar, unsafeLogVar, unsafeCov, unsafeLogCov, outboundType, isApplicable, dims, naturalParams, standardDistribution
using ForneyLab: SPLogNormalOutNPP, VBLogNormalOut

@testset "dims" begin
    @test dims(Distribution(Univariate, LogNormal, m=1.0, s=1.0)) == ()
end

@testset "vague" begin
    @test vague(LogNormal) == Distribution(Univariate, LogNormal, m=1.0, s=huge)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == exp(1.0 + 0.5*2.0)
    @test unsafeLogMean(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == 1.0
    @test unsafeVar(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == exp(2.0*1.0 + 2.0)*(exp(2.0) - 1.0)
    @test unsafeLogVar(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == 2.0
    @test unsafeCov(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == exp(2.0*1.0 + 2.0)*(exp(2.0) - 1.0)
    @test unsafeLogCov(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == 2.0
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, LogNormal, m=1.2, s=0.5),2), -1.522411904058978)
end

@testset "Gamma approximatons to LogNormal" begin
    @test ForneyLab.laplace(Gamma, Distribution(Univariate, LogNormal, m=0.0, s=2.0)) == Distribution(Univariate, Gamma, a=0.5, b=0.5)
end

@testset "prod!" begin
    @test Distribution(Univariate, LogNormal, m=1.0, s=2.0) * Distribution(Univariate, PointMass, m=1.0) == Distribution(Univariate, PointMass, m=1.0)
    @test_throws Exception Distribution(Univariate, LogNormal, m=1.0, s=2.0) * Distribution(Univariate, PointMass, m=-1.0)
    @test Distribution(Univariate, LogNormal, m=0.0, s=2.0) * Distribution(Univariate, Gamma, a=3.0, b=4.0) == Distribution(Univariate, Gamma, a=2.5, b=4.5)
end

@testset "natural parameters" begin
    d = Distribution(Univariate, LogNormal, m=1.0, s=2.0)
    η = naturalParams(d)
    s = standardDistribution(Univariate, LogNormal, η=η)
    @test d.params[:m] == s.params[:m] # Test conversion consistency
    @test d.params[:s] == s.params[:s]

    x = [0.2, 1.0, 8.0]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, LogNormal, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
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
    @test isApplicable(VBLogNormalOut, [Nothing, Distribution, Distribution])
    @test !isApplicable(VBLogNormalOut, [Distribution, Distribution, Nothing])

    @test ruleVBLogNormalOut(nothing, Distribution(Univariate, PointMass, m=1.5), Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, LogNormal, m=1.5, s=3.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, LogNormal, m=1.0, s=2.0)) == averageEnergy(LogNormal, Distribution(Univariate, LogNormal, m=1.0, s=2.0), Distribution(Univariate, PointMass, m=1.0), Distribution(Univariate, PointMass, m=2.0))
end

end #module
