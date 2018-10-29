module GaussianTest

using Test
using ForneyLab
import LinearAlgebra: Diagonal

@testset "sample" begin
    # Univariate
    @test isa(sample(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0)), Float64)

    # Multivariate
    @test isa(sample(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.2, 2.7], v=[2.0 -0.5; -0.5 1.5])), Vector{Float64})
    @test isa(sample(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.2, 2.7], v=Diagonal([2.0, 1.5]))), Vector{Float64})
end

@testset "prod!" begin
    # Univariate
    @test ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) * ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=2.0)
    @test ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) * ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=2.0)
    @test ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) * ProbabilityDistribution(Univariate, PointMass, m=1.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)
    @test ProbabilityDistribution(Univariate, PointMass, m=1.0) * ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == ProbabilityDistribution(Univariate, PointMass, m=1.0)

    # Multivariate
    @test ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=diageye(2)) * ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=diageye(2)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=2.0*diageye(2))
    @test ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)) * ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=2.0*diageye(2))
    @test ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)) * ProbabilityDistribution(Multivariate, PointMass, m=ones(2)) == ProbabilityDistribution(Multivariate, PointMass, m=ones(2))
    @test ProbabilityDistribution(Multivariate, PointMass, m=ones(2)) * ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)) == ProbabilityDistribution(Multivariate, PointMass, m=ones(2))
end

end #module