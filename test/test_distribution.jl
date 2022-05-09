module DistributionTest

using Test
using ForneyLab
using ForneyLab: isProper, gaussianQuadrature, matches
using LinearAlgebra: Diagonal

@testset "<<" begin
    @test <<(Distribution{Univariate, Gaussian}, Distribution)
    @test <<(Distribution{Univariate, Gaussian{Moments}}, Distribution)
    @test <<(Distribution{Multivariate, Gaussian{Moments}}, Distribution)
    @test !<<(Nothing, Distribution)
end

@testset "sample" begin
    dist = Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    samples = sample(dist, 10)
    @test length(samples) == 10
end

@testset "Distribution constructor aliases" begin
    @test P(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test ProbabilityDistribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
end

@testset "Univariate" begin
    # Distribution should be parameterized on a node type (distribution family)
    gaussian = Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test isa(gaussian, Distribution{Univariate})
    @test isa(gaussian, Distribution{Univariate, Gaussian{Moments}})
    @test !isa(gaussian, Distribution{Multivariate})
    @test !isa(gaussian, Distribution{MatrixVariate})
    @test gaussian.params[:m] == 0.0
    @test gaussian.params[:v] == 1.0

    # PointMass should be defined as a special family
    point_mass = Distribution(Univariate, PointMass, m=0.0)
    @test isa(point_mass, Distribution{Univariate, PointMass})
    @test point_mass.params == Dict(:m=>0.0)
    @test isProper(point_mass)
    @test mean(point_mass) == 0.0
    @test var(point_mass) == 0.0

    # Remaining DeltaFactors should not be allowed as a distribution family
    @test_throws Exception Distribution(Univariate, Equality)
    @test_throws Exception Distribution(Univariate, Clamp, m=0.0)

    @test Distribution(Univariate, PointMass, m=0.0) != Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test Distribution(Univariate, PointMass, m=0.0) != Distribution(Univariate, PointMass, m=1.0)
end

@testset "Multivariate" begin
    gaussian = Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0))
    @test isa(gaussian, Distribution{Multivariate})
    @test isa(gaussian, Distribution{Multivariate, Gaussian{Moments}})
    @test !isa(gaussian, Distribution{Univariate})
    @test !isa(gaussian, Distribution{MatrixVariate})
    @test gaussian.params[:m] == [0.0]
    @test gaussian.params[:v] == mat(1.0)

    # PointMass should be defined as a special family
    point_mass = Distribution(Multivariate, PointMass, m=[0.0])
    @test isa(point_mass, Distribution{Multivariate, PointMass})
    @test point_mass.params == Dict(:m=>[0.0])
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0]
    @test var(point_mass) == [0.0]
    @test cov(point_mass) == transpose([0.0])
end

@testset "MatrixVariate" begin
    point_mass = Distribution(MatrixVariate, PointMass, m=transpose([0.0]))
    @test isa(point_mass, Distribution{MatrixVariate, PointMass})
    @test point_mass.params == Dict(:m=>transpose([0.0]))
    @test isProper(point_mass)
    @test mean(point_mass) == transpose([0.0])
end

@testset "convert" begin
    @test convert(Distribution{Multivariate, PointMass}, Distribution(Univariate, PointMass, m=1.0)) == Distribution(Multivariate, PointMass, m=[1.0])
    @test convert(Distribution{MatrixVariate, PointMass}, Distribution(Univariate, PointMass, m=1.0)) == Distribution(MatrixVariate, PointMass, m=mat(1.0))
    @test convert(Distribution{MatrixVariate, PointMass}, Distribution(Multivariate, PointMass, m=[1.0, 2.0])) == Distribution(MatrixVariate, PointMass, m=reshape([1.0, 2.0], 2, 1))
end

@testset "PointMass Distribution and Message construction" begin
    @test Distribution(Univariate, PointMass, m=0.2) == Distribution{Univariate, PointMass}(Dict(:m=>0.2))
    @test_throws Exception Distribution(Multivariate, PointMass, m=0.2)
    @test Distribution(Multivariate, PointMass, m=[0.2]) == Distribution{Multivariate, PointMass}(Dict(:m=>[0.2]))
    @test Distribution(MatrixVariate, PointMass, m=transpose([0.2])) == Distribution{MatrixVariate, PointMass}(Dict(:m=>transpose([0.2])))
    @test Distribution(PointMass, m=0.2) == Distribution{Univariate, PointMass}(Dict(:m=>0.2))
    @test Distribution(Univariate, PointMass) == Distribution{Univariate, PointMass}(Dict(:m=>1.0))
    @test Distribution(PointMass) == Distribution{Univariate, PointMass}(Dict(:m=>1.0))
    @test Distribution(Multivariate, PointMass) == Distribution{Multivariate, PointMass}(Dict(:m=>[1.0]))
    @test Distribution(MatrixVariate, PointMass) == Distribution{MatrixVariate, PointMass}(Dict(:m=>mat(1.0)))
    @test Message(PointMass) == Message{PointMass, Univariate}(Distribution{Univariate, PointMass}(Dict(:m=>1.0)))
    @test Message(Univariate, PointMass) == Message{PointMass, Univariate}(Distribution{Univariate, PointMass}(Dict(:m=>1.0)))
    @test_throws Exception Message(Multivariate, PointMass, m=0.2)
end

@testset "Distribution and Message construction with functions" begin
    f(x) = x
    @test isa(Distribution(Univariate, Function, f=f).params[:f], Function)
    @test isa(Distribution(Univariate, Function, f=f), Distribution{Univariate, Function})
    @test isa(Distribution(Multivariate, Function, f=f), Distribution{Multivariate, Function})
    @test isa(Distribution(MatrixVariate, Function, f=f), Distribution{MatrixVariate, Function})
    @test isa(Distribution(Function, f=f), Distribution{Univariate, Function})
    @test isempty(vague(Function).params)
    @test isa(Message(Univariate, Function, f=()->()).dist.params[:f], Function)
    @test isa(Message(Univariate, Function), Message{Function, Univariate})
    @test isa(Message(Function), Message{Function, Univariate})
end

@testset "dims" begin
    @test dims(Distribution(Univariate, PointMass, m=0.0)) == ()
    @test dims(Distribution(Multivariate, PointMass, m=ones(2))) == (2,)
    @test dims(Distribution(MatrixVariate, PointMass, m=eye(2))) == (2,2)
end

@testset "gaussianQuadrature" begin
    @test gaussianQuadrature(x -> 1.0, m=-1.0, v=2.0) == 0.9999996733487053
    @test gaussianQuadrature(x -> x, m=-1.0, v=2.0) == -0.9999996733487055
    @test gaussianQuadrature(x -> x^2, m=-1.0, v=2.0) == 3.0000013419448215
end

end #module
