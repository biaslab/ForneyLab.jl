module SampleListTest

using Test
using Random
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, unsafeLogMean, unsafeMeanCov, unsafeMirroredLogMean, dims
using StatsFuns: betainvcdf
using SpecialFunctions: digamma

Random.seed!(1234)

@testset "SampleList ProbabilityDistribution construction" begin
    f_dummy(x) = x
    s = randn(10)
    @test dims(ProbabilityDistribution(Univariate, SampleList, s=s)) == 1
end

@testset "unsafeMean and unsafeVar" begin
    f_dummy(x) = x
    sigmoid(x) = 1/(1+exp(-x))
    
    s = randn(1000)
    @test abs(unsafeMean(ProbabilityDistribution(Univariate, SampleList, s=s, w=ones(1000)/1000)) - 0.0) < 0.1
    @test abs(unsafeVar(ProbabilityDistribution(Univariate, SampleList, s=s, w=ones(1000)/1000)) - 1.0) < 0.1
    
    sample_beta = x -> betainvcdf(0.5, 0.5, x)
    samples = sample_beta.(rand(1000))
    log_mean = digamma(0.5) - digamma(1.0)
    mirrored_log_mean = digamma(0.5) - digamma(1.0)
    
    @test abs(unsafeLogMean(ProbabilityDistribution(Univariate, SampleList, s=samples, w=ones(1000)/1000)) - log_mean) < 0.1
    @test abs(unsafeMirroredLogMean(ProbabilityDistribution(Univariate, SampleList, s=samples, w=ones(1000)/1000)) - mirrored_log_mean) < 0.1
end

@testset "prod!" begin
    s = randn(10000)
    p_dist = ProbabilityDistribution(Univariate, SampleList, s=s, w=ones(10000)/10000) * ProbabilityDistribution(GaussianMeanVariance, m=0.0, v=1.0)
    @test abs(unsafeMean(p_dist)) < 0.1
    @test abs(unsafeVar(p_dist) - 0.5) < 0.1
end

end
