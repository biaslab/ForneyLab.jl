module GaussianTest

using Test
using ForneyLab
using LinearAlgebra: Diagonal

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

#-------------
# Canonical Parameterization
#-------------

@testset "Exponential Family" begin
    # Univariate
    @test naturalParams(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.5, v=.5)) == naturalParams(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=5., w=2.))
    @test ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=5., w=2.) == standardDist(vague(GaussianMeanVariance), naturalParams(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.5, v=.5)))
    @test ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=5., w=2.) == standardMessage(vague(GaussianMeanVariance), naturalParams(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.5, v=.5))).dist
    @test isapprox(logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.5, v=.5), 1.0), logPdf(vague(GaussianMeanVariance), naturalParams(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.5, v=.5)), 1.0))
    @test isapprox(logPdf(ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=5., w=2.), 1.0), logPdf(vague(GaussianMeanVariance), naturalParams(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.5, v=.5)), 1.0))

    # Multivariate
    m_use, v_use = [-1.2, 2.3], [1.5 0.6;0.6 1.2]
    w_use = cholinv(v_use)
    xi_use = w_use*m_use
    @test naturalParams(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use)) == naturalParams(ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=xi_use, w=w_use))
    @test isapprox(xi_use, standardDist(vague(GaussianMeanVariance,2),naturalParams(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use))).params[:xi])
    @test isapprox(w_use, standardDist(vague(GaussianMeanVariance,2),naturalParams(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use))).params[:w])
    @test isapprox(xi_use, standardMessage(vague(GaussianMeanVariance,2),naturalParams(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use))).dist.params[:xi])
    @test isapprox(w_use, standardMessage(vague(GaussianMeanVariance,2),naturalParams(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use))).dist.params[:w])
    @test isapprox(logPdf(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use),[0.,1.]), logPdf(vague(GaussianMeanVariance,2),naturalParams(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_use, v=v_use)),[0.,1.]))
end

end #module