module ProbabilityDistributionTest

using Base.Test
import ForneyLab: ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, Gaussian, PointMass, Equality, mean, var, isProper, isValid, invalidate!, gaussianQuadrature

@testset "Univariate" begin
    # ProbabilityDistribution should be parameterized on a node type (distribution family)
    gaussian = Univariate(Gaussian, m=0.0, v=1.0)
    @test isa(gaussian, ProbabilityDistribution{Gaussian})
    @test isa(gaussian, Univariate{Gaussian})
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
    @test Univariate(PointMass, m=0.0) != Univariate(Gaussian, m=0.0, v=1.0)
    @test Univariate(PointMass, m=0.0) != Univariate(PointMass, m=1.0)
end

@testset "Multivariate" begin
    gaussian = Multivariate(Gaussian, m=[0.0], v=[1.0].')
    @test isa(gaussian, ProbabilityDistribution{Gaussian})
    @test isa(gaussian, Multivariate{Gaussian, 1})
    @test !isa(gaussian, Univariate)
    @test !isa(gaussian, MatrixVariate)
    @test gaussian.params == Dict(:m=>[0.0], :v=>[1.0].')

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
    point_mass = MatrixVariate(PointMass, m=[0.0].')
    @test isa(point_mass, MatrixVariate{PointMass, 1, 1})
    @test point_mass.params == Dict(:m=>[0.0].')
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0].'
end

@testset "isValid" begin
    # should check validity of distribution parameters
    @test isValid(Univariate(PointMass, m=1.0), :m)
    @test !isValid(Univariate(PointMass, m=NaN), :m)
    @test isValid(Multivariate(PointMass, m=[1.0, 2.0]), :m)
    @test !isValid(Multivariate(PointMass, m=[NaN, 2.0]), :m)
    @test isValid(MatrixVariate(PointMass, m=eye(2)), :m)
    @test !isValid(MatrixVariate(PointMass, m=[NaN 1.0; 2.0 3.0]), :m)
    @test isValid(MatrixVariate(PointMass, m=Diagonal([1.0, 2.0])), :m)
    @test !isValid(MatrixVariate(PointMass, m=Diagonal([NaN, 2.0])), :m)
end

@testset "invalidate!" begin
    # should invalidate vectors and matrices
    d = Multivariate(PointMass, m=[1.0, 2.0])
    @test !isValid(invalidate!(d, :m), :m)
    d = MatrixVariate(PointMass, m=[NaN 1.0; 2.0 3.0])
    @test !isValid(invalidate!(d, :m), :m)
    d = MatrixVariate(PointMass, m=Diagonal([1.0, 2.0]))
    @test !isValid(invalidate!(d, :m), :m)
end

@testset "gaussianQuadrature" begin
    @test gaussianQuadrature(x -> 1.0, m=-1.0, v=2.0) == 0.9999996733487053
    @test gaussianQuadrature(x -> x, m=-1.0, v=2.0) == -0.9999996733487055
    @test gaussianQuadrature(x -> x^2, m=-1.0, v=2.0) == 3.0000013419448215
end

end #module