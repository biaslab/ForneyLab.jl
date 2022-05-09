module GaussianTest

using Test
using ForneyLab
using ForneyLab: naturalParams, standardDistribution
using LinearAlgebra: Diagonal

@testset "Gaussian distribution constructor aliases" begin
    @test Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0) == Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0)
    @test Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0)
end

@testset "sample" begin
    # Univariate
    @test isa(sample(Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0)), Float64)

    # Multivariate
    @test isa(sample(Distribution(Multivariate, Gaussian{Moments}, m=[1.2, 2.7], v=[2.0 -0.5; -0.5 1.5])), Vector{Float64})
    @test isa(sample(Distribution(Multivariate, Gaussian{Moments}, m=[1.2, 2.7], v=Diagonal([2.0, 1.5]))), Vector{Float64})
end

@testset "prod!" begin
    # Univariate
    @test Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0) * Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=2.0)
    @test Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) * Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=2.0)
    @test Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) * Distribution(Univariate, PointMass, m=1.0) == Distribution(Univariate, PointMass, m=1.0)
    @test Distribution(Univariate, PointMass, m=1.0) * Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Distribution(Univariate, PointMass, m=1.0)

    # Multivariate
    @test Distribution(Multivariate, Gaussian{Canonical}, xi=zeros(2), w=diageye(2)) * Distribution(Multivariate, Gaussian{Canonical}, xi=zeros(2), w=diageye(2)) == Distribution(Multivariate, Gaussian{Canonical}, xi=zeros(2), w=2.0*diageye(2))
    @test Distribution(Multivariate, Gaussian{Moments}, m=zeros(2), v=diageye(2)) * Distribution(Multivariate, Gaussian{Moments}, m=zeros(2), v=diageye(2)) == Distribution(Multivariate, Gaussian{Canonical}, xi=zeros(2), w=2.0*diageye(2))
    @test Distribution(Multivariate, Gaussian{Moments}, m=zeros(2), v=diageye(2)) * Distribution(Multivariate, PointMass, m=ones(2)) == Distribution(Multivariate, PointMass, m=ones(2))
    @test Distribution(Multivariate, PointMass, m=ones(2)) * Distribution(Multivariate, Gaussian{Moments}, m=zeros(2), v=diageye(2)) == Distribution(Multivariate, PointMass, m=ones(2))
end

@testset "natural parameters" begin
    # Univariate
    d = Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=5.0)
    η = naturalParams(d)
    s = standardDistribution(Univariate, Gaussian, η=η)
    @test d.params[:xi] == s.params[:xi] # Test conversion consistency
    @test d.params[:w] == s.params[:w]

    x = [-0.1, 2.0, 15.0]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, Gaussian, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency

    # Multivariate
    d = Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0, 5.0], w=[2.0 0.1; 0.1 3.0])
    η = naturalParams(d)
    s = standardDistribution(Multivariate, Gaussian, η=η)
    @test d.params[:xi] == s.params[:xi] # Test conversion consistency
    @test d.params[:w] == s.params[:w]

    x = [[-0.1, 2.0], [-3.0, 5.0]]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Multivariate, Gaussian, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
end

end #module