module GaussianWeightedMeanPrecisionTest

using Test
using ForneyLab
import ForneyLab: isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=ones(2), w=diageye(2))) == 2
end

@testset "vague" begin
    @test vague(GaussianWeightedMeanPrecision) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
    @test vague(GaussianWeightedMeanPrecision, 2) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=tiny*eye(2))
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

@testset "convert" begin
    @test convert(ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.5, w=4.0)) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)
    @test convert(ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.5, v=0.25)) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)
    @test convert(ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.5], w=mat(4.0))) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))
    @test convert(ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.5], v=mat(0.25))) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))
end

end # module