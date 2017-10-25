module ProbabilityDistributionTest

using Base.Test
import ForneyLab: ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, GaussianMeanVariance, PointMass, Equality, mean, var, isProper, gaussianQuadrature

@testset "Univariate" begin
    # ProbabilityDistribution should be parameterized on a node type (distribution family)
    gaussian = Univariate(GaussianMeanVariance, m=0.0, v=1.0)
    @test_throws Exception Univariate(GaussianMeanVariance, m=0.0, V=1.0)
    @test isa(gaussian, ProbabilityDistribution{GaussianMeanVariance})
    @test isa(gaussian, Univariate{GaussianMeanVariance})
    @test !isa(gaussian, Multivariate)
    @test !isa(gaussian, MatrixVariate)
    @test gaussian.params == Dict(:m=>0.0, :v=>1.0)

    # PointMass should be defined as a special family 
    point_mass = Univariate(PointMass, m=0.0)
    @test isa(point_mass, Univariate{PointMass})
    @test point_mass.params == Dict(:m=>0.0)
    @test isProper(point_mass)
    @test mean(point_mass) == 0.0
    @test var(point_mass) == 0.0

    # Remaining DeltaFactors should not be allowed as a distribution family
    @test_throws Exception Univariate(Equality)
    @test_throws Exception Univariate(Clamp, m=0.0)

    @test Univariate(PointMass, m=0.0) == Univariate(PointMass, m=0.0)
    @test Univariate(PointMass, m=0.0) != Univariate(GaussianMeanVariance, m=0.0, v=1.0)
    @test Univariate(PointMass, m=0.0) != Univariate(PointMass, m=1.0)
end

@testset "Multivariate" begin
    gaussian = Multivariate(GaussianMeanVariance, m=[0.0], V=[1.0].')
    @test_throws Exception Multivariate(GaussianMeanVariance, m=0.0, V=1.0)
    @test_throws Exception Multivariate(GaussianMeanVariance, m=[0.0], v=[1.0].')
    @test isa(gaussian, ProbabilityDistribution{GaussianMeanVariance})
    @test isa(gaussian, Multivariate{GaussianMeanVariance, 1})
    @test !isa(gaussian, Univariate)
    @test !isa(gaussian, MatrixVariate)
    @test gaussian.params == Dict(:m=>[0.0], :V=>[1.0].')

    # PointMass should be defined as a special family 
    point_mass = Multivariate(PointMass, m=[0.0])
    @test isa(point_mass, Multivariate{PointMass, 1})
    @test point_mass.params == Dict(:m=>[0.0])
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0]
    @test var(point_mass) == [0.0]
    @test cov(point_mass) == [0.0].'
end

@testset "MatrixVariate" begin
    wishart = MatrixVariate(Wishart, nu=1.0, V=[1.0].')
    @test isa(wishart, ProbabilityDistribution{Wishart})
    @test isa(wishart, MatrixVariate{Wishart, 1, 1})
    @test !isa(wishart, Univariate)
    @test !isa(wishart, Multivariate)
    @test wishart.params == Dict(:nu=>1.0, :V=>[1.0].')

    # PointMass should be defined as a special family 
    point_mass = MatrixVariate(PointMass, M=[0.0].')
    @test isa(point_mass, MatrixVariate{PointMass, 1, 1})
    @test point_mass.params == Dict(:M=>[0.0].')
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0].'
end

@testset "isValid" begin
    # should check validity of distribution parameters
    @test isValid(Univariate(PointMass, m=1.0), :m)
    @test !isValid(Univariate(PointMass, m=NaN), :m)
    @test isValid(Multivariate(PointMass, m=[1.0, 2.0]), :m)
    @test !isValid(Multivariate(PointMass, m=[NaN, 2.0]), :m)
    @test isValid(MatrixVariate(PointMass, M=eye(2)), :M)
    @test !isValid(MatrixVariate(PointMass, M=[NaN 1.0; 2.0 3.0]), :M)
    @test isValid(MatrixVariate(PointMass, M=Diagonal([1.0, 2.0])), :M)
    @test !isValid(MatrixVariate(PointMass, M=Diagonal([NaN, 2.0])), :M)
end

@testset "invalidate!" begin
    # should invalidate vectors and matrices
    A = [1.0, 2.0]
    @test isValid(invalidate!(A)) == false
    @test isValid(A) == false
    A = [1.0 2.0; 3.0 4.0]
    @test isValid(invalidate!(A)) == false
    @test isValid(A) == false
    A = Diagonal([1.0, 2.0])
    @test isValid(invalidate!(A)) == false
    @test isValid(A) == false
end

@testset "gaussianQuadrature" begin
    @test gaussianQuadrature(x -> 1.0, m=-1.0, v=2.0) == 0.9999996733487053
    @test gaussianQuadrature(x -> x, m=-1.0, v=2.0) == -0.9999996733487055
    @test gaussianQuadrature(x -> x^2, m=-1.0, v=2.0) == 3.0000013419448215
end

end #module