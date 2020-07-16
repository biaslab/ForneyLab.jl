module SampleListTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeCov, unsafeVar, dims
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

end
