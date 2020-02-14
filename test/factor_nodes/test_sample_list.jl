module SampleListTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, unsafeLogMean, unsafeMeanCov, unsafeMirroredLogMean, dims
using StatsFuns: betainvcdf
using SpecialFunctions: digamma

@testset "SampleList ProbabilityDistribution construction" begin
    @test_throws Exception ProbabilityDistribution(Multivariate, SampleList)

    f_dummy(x) = x
    s = randn(10)
    @test dims(ProbabilityDistribution(Univariate, SampleList, s=s)) == 1
end

@testset "unsafeMean and unsafeVar" begin
    f_dummy(x) = x
    sigmoid(x) = 1/(1+exp(-x))
    
    s = randn(1000)
    @test abs(unsafeMean(ProbabilityDistribution(Univariate, SampleList, s=s)) - 0.0) < 0.1
    @test abs(unsafeVar(ProbabilityDistribution(Univariate, SampleList, s=s)) - 1.0) < 0.1
    
    sample_beta = x -> betainvcdf(0.5, 0.5, x)
    samples = sample_beta.(rand(1000))
    log_mean = digamma(0.5) - digamma(1.0)
    mirrored_log_mean = digamma(0.5) - digamma(1.0)
    
    @test abs(unsafeLogMean(ProbabilityDistribution(Univariate, SampleList, s=samples)) - log_mean) < 0.1
    @test abs(unsafeMirroredLogMean(ProbabilityDistribution(Univariate, SampleList, s=samples)) - mirrored_log_mean) < 0.1
end

@testset "prod!" begin
    s = randn(10000)
    p_dist = ProbabilityDistribution(Univariate, SampleList, s=s) * ProbabilityDistribution(GaussianMeanVariance, m=0.0, v=1.0)
    @test abs(unsafeMean(p_dist)) < 0.1
    @test abs(unsafeVar(p_dist) - 0.5) < 0.1
end

end
