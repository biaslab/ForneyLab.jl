module SampleListTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeCov, unsafeVar, unsafeLogMean, unsafeMirroredLogMean, dims
using ForneyLab: SPSampleListOutNPP

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.5, 0.5])) == 1
    @test dims(ProbabilityDistribution(Multivariate, SampleList, s=[[0.0], [1.0]], w=[0.5, 0.5])) == 1
    @test dims(ProbabilityDistribution(MatrixVariate, SampleList, s=[mat(0.0), mat(1.0)], w=[0.5, 0.5])) == (1,1)
end

@testset "unsafeMean" begin
    @test unsafeMean(ProbabilityDistribution(Univariate, SampleList, s=[1.0, 1.0], w=[0.5, 0.5])) == 1.0
    @test unsafeMean(ProbabilityDistribution(Multivariate, SampleList, s=[[1.0], [1.0]], w=[0.5, 0.5])) == [1.0]
    @test unsafeMean(ProbabilityDistribution(MatrixVariate, SampleList, s=[mat(1.0), mat(1.0)], w=[0.5, 0.5])) == mat(1.0)
end

@testset "unsafeVar" begin
    @test unsafeVar(ProbabilityDistribution(Univariate, SampleList, s=[1.0, 1.0], w=[0.5, 0.5])) == 0.0
    @test unsafeVar(ProbabilityDistribution(Multivariate, SampleList, s=[[1.0], [1.0]], w=[0.5, 0.5])) == [0.0]
end

@testset "unsafeCov" begin
    @test unsafeCov(ProbabilityDistribution(Univariate, SampleList, s=[1.0, 1.0], w=[0.5, 0.5])) == 0.0
    @test unsafeCov(ProbabilityDistribution(Multivariate, SampleList, s=[[1.0], [1.0]], w=[0.5, 0.5])) == mat(0.0)
    @test unsafeCov(ProbabilityDistribution(MatrixVariate, SampleList, s=[eye(2), eye(2)], w=[0.5, 0.5])) == zeros(4,4)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.5, 0.5]) * ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == ProbabilityDistribution(Univariate, SampleList, s=[0.0, 1.0], w=[0.6224593312018546, 0.37754066879814546])
    @test ProbabilityDistribution(Multivariate, SampleList, s=[[0.0], [1.0]], w=[0.5, 0.5]) * ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)) == ProbabilityDistribution(Multivariate, SampleList, s=[[0.0], [1.0]], w=[0.6224593312018546, 0.37754066879814546])
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

#------------
# Integration
#------------
using StatsFuns: betainvcdf
using SpecialFunctions: digamma
using LinearAlgebra: norm
using Random

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

@testset "unsafeMean, unsafeVar,unsafeLogMean, unsafeMirroredLogMean numeric test" begin
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

end
