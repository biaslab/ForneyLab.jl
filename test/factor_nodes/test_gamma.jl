module GammaTest

using Test
using ForneyLab
using ForneyLab: prod!, unsafeMean, unsafeVar, outboundType, isApplicable, dims, naturalParams, standardDistribution
using ForneyLab: SPGammaOutNPP, VBGammaOut, VBGammaA, VBGammaB

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=1.0)) == ()
end

@testset "vague" begin
    @test vague(Gamma) == ProbabilityDistribution(Univariate, Gamma, a=1.0, b=tiny)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0) * ProbabilityDistribution(Univariate, Gamma, a=3.0, b=4.0) == ProbabilityDistribution(Univariate, Gamma, a=3.0, b=6.0)
    @test ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0) * ProbabilityDistribution(Univariate, PointMass, m=1.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)
    @test ProbabilityDistribution(Univariate, PointMass, m=1.0) * ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)
    @test_throws Exception ProbabilityDistribution(Univariate, PointMass, m=-1.0) * ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == 0.5
    @test unsafeVar(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == 0.25
end

@testset "natural parameters" begin
    d = ProbabilityDistribution(Univariate, Gamma, a=2.0, b=5.0)
    η = naturalParams(d)
    s = standardDistribution(Univariate, Gamma, η=η)
    @test d.params[:a] == s.params[:a] # Test conversion consistency
    @test d.params[:b] == s.params[:b]

    x = [0.1, 2.0, 15.0]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, Gamma, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
end


#-------------
# Update rules
#-------------

@testset "SPGammaOutNPP" begin
    @test SPGammaOutNPP <: SumProductRule{Gamma}
    @test outboundType(SPGammaOutNPP) == Message{Gamma}
    @test isApplicable(SPGammaOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPGammaOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gamma, a=1.0, b=2.0)
end

@testset "VBGammaOut" begin
    @test VBGammaOut <: NaiveVariationalRule{Gamma}
    @test outboundType(VBGammaOut) == Message{Gamma}
    @test isApplicable(VBGammaOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBGammaOut, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBGammaOut(nothing, ProbabilityDistribution(Univariate, PointMass, m=1.5), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gamma, a=1.5, b=3.0)
end

@testset "VBGammaA" begin
    @test VBGammaA <: NaiveVariationalRule{Gamma}
    @test outboundType(VBGammaA) == Message{Function}
    @test !isApplicable(VBGammaA, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test isApplicable(VBGammaA, [ProbabilityDistribution, Nothing, ProbabilityDistribution])
end

@testset "VBGammaB" begin
    @test VBGammaB <: NaiveVariationalRule{Gamma}
    @test outboundType(VBGammaB) == Message{Gamma}
    @test !isApplicable(VBGammaB, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test isApplicable(VBGammaB, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBGammaB(ProbabilityDistribution(Univariate, PointMass, m=1.5), ProbabilityDistribution(Univariate, PointMass, m=3.0), nothing) == Message(Univariate, Gamma, a=4.0, b=1.5)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == averageEnergy(Gamma, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, PointMass, m=1.0), ProbabilityDistribution(Univariate, PointMass, m=2.0))
end

end #module
