module ProbabilityDistributionTest

using Base.Test
import ForneyLab: ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, Gaussian, PointMass, Equality, mean, var, isProper, isValid, invalidate!, gaussianQuadrature, dims, diageye, matches

@testset "matches" begin
    @test matches(ProbabilityDistribution{Univariate, Gaussian}, ProbabilityDistribution{Univariate, Gaussian})
    @test !matches(ProbabilityDistribution{Univariate, Gaussian}, ProbabilityDistribution{Univariate, PointMass})
    @test matches(ProbabilityDistribution{Univariate, Gaussian}, ProbabilityDistribution{Univariate})
    @test matches(ProbabilityDistribution{Univariate, Gaussian}, ProbabilityDistribution)
    @test matches(ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate})
    @test !matches(ProbabilityDistribution{Univariate}, ProbabilityDistribution{Multivariate})
end

@testset "Univariate" begin
    # ProbabilityDistribution should be parameterized on a node type (distribution family)
    gaussian = Univariate(Gaussian, m=0.0, v=1.0)
    @test isa(gaussian, ProbabilityDistribution{Univariate})
    @test isa(gaussian, ProbabilityDistribution{Univariate, Gaussian})
    @test !isa(gaussian, ProbabilityDistribution{Multivariate})
    @test !isa(gaussian, ProbabilityDistribution{MatrixVariate})
    @test gaussian.params == Dict(:m=>0.0, :v=>1.0)

    # Test ProbabilityDistribution constructor
    @test ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0) == gaussian 

    # PointMass should be defined as a special family 
    point_mass = Univariate(PointMass, m=0.0)
    @test isa(point_mass, ProbabilityDistribution{Univariate, PointMass})
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
    @test isa(gaussian, ProbabilityDistribution{Multivariate})
    @test isa(gaussian, ProbabilityDistribution{Multivariate, Gaussian})
    @test !isa(gaussian, Univariate)
    @test !isa(gaussian, MatrixVariate)
    @test gaussian.params == Dict(:m=>[0.0], :v=>[1.0].')

    # PointMass should be defined as a special family 
    point_mass = Multivariate(PointMass, m=[0.0])
    @test isa(point_mass, ProbabilityDistribution{Multivariate, PointMass})
    @test point_mass.params == Dict(:m=>[0.0])
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0]
    @test var(point_mass) == [0.0]
    @test cov(point_mass) == [0.0].'
end

@testset "MatrixVariate" begin
    point_mass = MatrixVariate(PointMass, m=[0.0].')
    @test isa(point_mass, ProbabilityDistribution{MatrixVariate, PointMass})
    @test point_mass.params == Dict(:m=>[0.0].')
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0].'
end

@testset "dims" begin
    @test dims(Univariate(PointMass, m=0.0)) == 1
    @test dims(Multivariate(PointMass, m=ones(2))) == 2
    @test dims(MatrixVariate(PointMass, m=eye(2))) == (2, 2)
    @test dims(MatrixVariate(PointMass, m=diageye(2))) == (2, 2)
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