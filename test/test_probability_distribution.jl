module ProbabilityDistributionTest

using Base.Test
import ForneyLab: ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, Gaussian, PointMass, Equality, mean, var, mat, isProper, isValid, invalidate!, gaussianQuadrature, dims, diageye, matches, Message

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
    gaussian = ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)
    @test isa(gaussian, ProbabilityDistribution{Univariate})
    @test isa(gaussian, ProbabilityDistribution{Univariate, Gaussian})
    @test !isa(gaussian, ProbabilityDistribution{Multivariate})
    @test !isa(gaussian, ProbabilityDistribution{MatrixVariate})
    @test gaussian.params[:m] == 0.0
    @test gaussian.params[:v] == 1.0
    @test isnan(gaussian.params[:w])
    @test isnan(gaussian.params[:xi])

    # PointMass should be defined as a special family 
    point_mass = ProbabilityDistribution(Univariate, PointMass, m=0.0)
    @test isa(point_mass, ProbabilityDistribution{Univariate, PointMass})
    @test point_mass.params == Dict(:m=>0.0)
    @test isProper(point_mass)
    @test mean(point_mass) == 0.0
    @test var(point_mass) == 0.0

    # Remaining DeltaFactors should not be allowed as a distribution family
    @test_throws Exception ProbabilityDistribution(Univariate, Equality)
    @test_throws Exception ProbabilityDistribution(Univariate, Clamp, m=0.0)

    @test ProbabilityDistribution(Univariate, PointMass, m=0.0) != ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)
    @test ProbabilityDistribution(Univariate, PointMass, m=0.0) != ProbabilityDistribution(Univariate, PointMass, m=1.0)
end

@testset "Multivariate" begin
    gaussian = ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=mat(1.0))
    @test isa(gaussian, ProbabilityDistribution{Multivariate})
    @test isa(gaussian, ProbabilityDistribution{Multivariate, Gaussian})
    @test !isa(gaussian, ProbabilityDistribution{Univariate})
    @test !isa(gaussian, ProbabilityDistribution{MatrixVariate})
    @test gaussian.params[:m] == [0.0]
    @test gaussian.params[:v] == mat(1.0)
    @test isnan(gaussian.params[:w][1])
    @test isnan(gaussian.params[:xi][1])

    # PointMass should be defined as a special family 
    point_mass = ProbabilityDistribution(Multivariate, PointMass, m=[0.0])
    @test isa(point_mass, ProbabilityDistribution{Multivariate, PointMass})
    @test point_mass.params == Dict(:m=>[0.0])
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0]
    @test var(point_mass) == [0.0]
    @test cov(point_mass) == [0.0].'
end

@testset "MatrixVariate" begin
    point_mass = ProbabilityDistribution(MatrixVariate, PointMass, m=[0.0].')
    @test isa(point_mass, ProbabilityDistribution{MatrixVariate, PointMass})
    @test point_mass.params == Dict(:m=>[0.0].')
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0].'
end

@testset "PointMass ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Univariate, PointMass, m=0.2) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>0.2))
    @test_throws Exception ProbabilityDistribution(Multivariate, PointMass, m=0.2)
    @test ProbabilityDistribution(Multivariate, PointMass, m=[0.2]) == ProbabilityDistribution{Multivariate, PointMass}(Dict(:m=>[0.2]))
    @test ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2].') == ProbabilityDistribution{MatrixVariate, PointMass}(Dict(:m=>[0.2].'))
    @test ProbabilityDistribution(PointMass, m=0.2) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>0.2))
    @test ProbabilityDistribution(Univariate, PointMass) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0))
    @test ProbabilityDistribution(PointMass) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0))
    @test ProbabilityDistribution(Multivariate, PointMass) == ProbabilityDistribution{Multivariate, PointMass}(Dict(:m=>[1.0]))
    @test ProbabilityDistribution(MatrixVariate, PointMass) == ProbabilityDistribution{MatrixVariate, PointMass}(Dict(:m=>mat(1.0)))
    @test Message(PointMass) == Message{PointMass, Univariate}(ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0)))
    @test Message(Univariate, PointMass) == Message{PointMass, Univariate}(ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0)))
    @test_throws Exception Message(Multivariate, PointMass, m=0.2)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, PointMass, m=0.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, PointMass, m=ones(2))) == 2
    @test dims(ProbabilityDistribution(MatrixVariate, PointMass, m=eye(2))) == (2, 2)
    @test dims(ProbabilityDistribution(MatrixVariate, PointMass, m=diageye(2))) == (2, 2)
end

@testset "isValid" begin
    # should check validity of distribution parameters
    @test isValid(ProbabilityDistribution(Univariate, PointMass, m=1.0), :m)
    @test !isValid(ProbabilityDistribution(Univariate, PointMass, m=NaN), :m)
    @test isValid(ProbabilityDistribution(Multivariate, PointMass, m=[1.0, 2.0]), :m)
    @test !isValid(ProbabilityDistribution(Multivariate, PointMass, m=[NaN, 2.0]), :m)
    @test isValid(ProbabilityDistribution(MatrixVariate, PointMass, m=eye(2)), :m)
    @test !isValid(ProbabilityDistribution(MatrixVariate, PointMass, m=[NaN 1.0; 2.0 3.0]), :m)
    @test isValid(ProbabilityDistribution(MatrixVariate, PointMass, m=Diagonal([1.0, 2.0])), :m)
    @test !isValid(ProbabilityDistribution(MatrixVariate, PointMass, m=Diagonal([NaN, 2.0])), :m)
end

@testset "invalidate!" begin
    # should invalidate vectors and matrices
    d = ProbabilityDistribution(Multivariate, PointMass, m=[1.0, 2.0])
    @test !isValid(invalidate!(d, :m), :m)
    d = ProbabilityDistribution(MatrixVariate, PointMass, m=[NaN 1.0; 2.0 3.0])
    @test !isValid(invalidate!(d, :m), :m)
    d = ProbabilityDistribution(MatrixVariate, PointMass, m=Diagonal([1.0, 2.0]))
    @test !isValid(invalidate!(d, :m), :m)
end

@testset "gaussianQuadrature" begin
    @test gaussianQuadrature(x -> 1.0, m=-1.0, v=2.0) == 0.9999996733487053
    @test gaussianQuadrature(x -> x, m=-1.0, v=2.0) == -0.9999996733487055
    @test gaussianQuadrature(x -> x^2, m=-1.0, v=2.0) == 3.0000013419448215
end

end #module