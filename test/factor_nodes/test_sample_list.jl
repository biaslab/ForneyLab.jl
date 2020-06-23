module SampleListTest

using Test
using Random
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeCov, unsafeVar, unsafeLogMean, unsafeMeanCov, unsafeMirroredLogMean, dims
using ForneyLab: SPSampleListOutNPP
using StatsFuns: betainvcdf
using SpecialFunctions: digamma
using LinearAlgebra: norm

Random.seed!(1234)

@testset "SampleList ProbabilityDistribution construction" begin
    f_dummy(x) = x
    s = randn(10)
    @test dims(ProbabilityDistribution(Univariate, SampleList, s=s)) == 1

    m = Vector{Matrix}(undef,10)
    for i=1:10
        m[i] = randn(3,4)
    end
    @test dims(ProbabilityDistribution(MatrixVariate,SampleList,s=m)) == (3,4)
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

    m = Vector{Matrix}(undef,100000)
    mt = Vector{Matrix}(undef,100000)
    for i=1:100000
        m[i] = randn(3,4)
        mt[i] = transpose(m[i])
    end
    @test norm(unsafeMean(ProbabilityDistribution(MatrixVariate,SampleList,s=m,w=ones(100000)/100000)) - mean(m)) < 0.3
    @test norm(unsafeCov(ProbabilityDistribution(MatrixVariate,SampleList,s=m,w=ones(100000)/100000)) - kron(cov(m),cov(mt))) < 0.3
end

@testset "prod!" begin
    s = randn(10000)
    p_dist = ProbabilityDistribution(Univariate, SampleList, s=s, w=ones(10000)/10000) * ProbabilityDistribution(GaussianMeanVariance, m=0.0, v=1.0)
    @test abs(unsafeMean(p_dist)) < 0.1
    @test abs(unsafeVar(p_dist) - 0.5) < 0.1
end


#-------------
# Update rules
#-------------

@testset "SPSampleListOutNPP" begin
    @test SPSampleListOutNPP <: SumProductRule{SampleList}
    @test outboundType(SPSampleListOutNPP) == Message{SampleList}
    @test isApplicable(SPSampleListOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPSampleListOutNPP(nothing, Message(Multivariate, PointMass, m=[0.0,2.0,5.2]), Message(Multivariate, PointMass, m=[0.4,0.5,0.1])) == Message(Univariate, SampleList, s=[0.0,2.0,5.2], w=[0.4,0.5,0.1])
    @test ruleSPSampleListOutNPP(nothing, Message(Multivariate, PointMass, m=[[0.0,2.0],[5.2,-0.6]]), Message(Multivariate, PointMass, m=[0.4,0.6])) == Message(Multivariate, SampleList, s=[[0.0,2.0],[5.2,-0.6]], w=[0.4,0.6])
end

end
