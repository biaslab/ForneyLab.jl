module GaussianWeightedMeanPrecisionTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianWeightedMeanPrecisionOutNPP, VBGaussianWeightedMeanPrecisionOut

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=ones(2), w=diageye(2))) == 2
end

@testset "vague" begin
    @test vague(GaussianWeightedMeanPrecision) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
    @test vague(GaussianWeightedMeanPrecision, 2) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=tiny*eye(2))
    @test vague(GaussianWeightedMeanPrecision, (2,)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=tiny*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0))
    @test !isProper(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=-1.0))

    # Multivariate
    @test isProper(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0)))
    @test isProper(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=ones(2), w=diageye(2)))
    @test !isProper(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)
    @test ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)

    # Multivariate
    @test ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0))
    @test ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.5
    @test unsafeMode(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.5
    @test unsafeVar(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.25
    @test unsafeCov(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.25
    @test unsafeMeanCov(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == (0.5, 0.25)
    @test unsafePrecision(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 4.0
    @test unsafeWeightedMean(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 2.0
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == (2.0, 4.0)

    # Multivariate
    @test unsafeMean(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [0.5]
    @test unsafeMode(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [0.5]
    @test unsafeVar(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [0.25]
    @test unsafeCov(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == mat(0.25)
    @test unsafeMeanCov(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == ([0.5], mat(0.25))
    @test unsafePrecision(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == mat(4.0)
    @test unsafeWeightedMean(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == ([2.0], mat(4.0))
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=2.0), 1.0), -0.5723649429247)
    @test isapprox(logPdf(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0, 2.0], w=[2.0 0.0; 0.0 2.0]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.5, w=4.0)) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)
    @test convert(ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.5, v=0.25)) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)
    @test convert(ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.5], w=mat(4.0))) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))
    @test convert(ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.5], v=mat(0.25))) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))
    @test convert(ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0))
end

#-------------
# Update rules
#-------------

@testset "SPGaussianWeightedMeanPrecisionOutNPP" begin
    @test SPGaussianWeightedMeanPrecisionOutNPP <: SumProductRule{GaussianWeightedMeanPrecision}
    @test outboundType(SPGaussianWeightedMeanPrecisionOutNPP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPGaussianWeightedMeanPrecisionOutNPP, [Nothing, Message{PointMass}, Message{PointMass}]) 
    @test !isApplicable(SPGaussianWeightedMeanPrecisionOutNPP, [Message{PointMass}, Nothing, Message{PointMass}]) 

    @test ruleSPGaussianWeightedMeanPrecisionOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0)
    @test ruleSPGaussianWeightedMeanPrecisionOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0))
end

@testset "VBGaussianWeightedMeanPrecisionOut" begin
    @test VBGaussianWeightedMeanPrecisionOut <: NaiveVariationalRule{GaussianWeightedMeanPrecision}
    @test outboundType(VBGaussianWeightedMeanPrecisionOut) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(VBGaussianWeightedMeanPrecisionOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBGaussianWeightedMeanPrecisionOut, [ProbabilityDistribution, ProbabilityDistribution, Nothing]) 

    @test ruleVBGaussianWeightedMeanPrecisionOut(nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=3.0)
    @test ruleVBGaussianWeightedMeanPrecisionOut(nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(3.0))
end

end # module