module ProbabilityDistributionTest

using Test
import ForneyLab: ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, Gaussian, PointMass, Equality, mean, var, mat, isProper, isValid, invalidate!, gaussianQuadrature, dims, eye, diageye, matches, Message
import LinearAlgebra: Diagonal
using ForneyLab

@testset "matches" begin
    @test matches(ProbabilityDistribution{Univariate, Gaussian}, ProbabilityDistribution)
    @test matches(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution)
    @test matches(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, ProbabilityDistribution)
    @test !matches(Nothing, ProbabilityDistribution)
end

@testset "Univariate" begin
    # ProbabilityDistribution should be parameterized on a node type (distribution family)
    gaussian = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    @test isa(gaussian, ProbabilityDistribution{Univariate})
    @test isa(gaussian, ProbabilityDistribution{Univariate, GaussianMeanVariance})
    @test !isa(gaussian, ProbabilityDistribution{Multivariate})
    @test !isa(gaussian, ProbabilityDistribution{MatrixVariate})
    @test gaussian.params[:m] == 0.0
    @test gaussian.params[:v] == 1.0

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

    @test ProbabilityDistribution(Univariate, PointMass, m=0.0) != ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    @test ProbabilityDistribution(Univariate, PointMass, m=0.0) != ProbabilityDistribution(Univariate, PointMass, m=1.0)
end

@testset "Multivariate" begin
    gaussian = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    @test isa(gaussian, ProbabilityDistribution{Multivariate})
    @test isa(gaussian, ProbabilityDistribution{Multivariate, GaussianMeanVariance})
    @test !isa(gaussian, ProbabilityDistribution{Univariate})
    @test !isa(gaussian, ProbabilityDistribution{MatrixVariate})
    @test gaussian.params[:m] == [0.0]
    @test gaussian.params[:v] == mat(1.0)

    # PointMass should be defined as a special family
    point_mass = ProbabilityDistribution(Multivariate, PointMass, m=[0.0])
    @test isa(point_mass, ProbabilityDistribution{Multivariate, PointMass})
    @test point_mass.params == Dict(:m=>[0.0])
    @test isProper(point_mass)
    @test mean(point_mass) == [0.0]
    @test var(point_mass) == [0.0]
    @test cov(point_mass) == transpose([0.0])
end

@testset "MatrixVariate" begin
    point_mass = ProbabilityDistribution(MatrixVariate, PointMass, m=transpose([0.0]))
    @test isa(point_mass, ProbabilityDistribution{MatrixVariate, PointMass})
    @test point_mass.params == Dict(:m=>transpose([0.0]))
    @test isProper(point_mass)
    @test mean(point_mass) == transpose([0.0])
end

@testset "PointMass ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Univariate, PointMass, m=0.2) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>0.2))
    @test_throws Exception ProbabilityDistribution(Multivariate, PointMass, m=0.2)
    @test ProbabilityDistribution(Multivariate, PointMass, m=[0.2]) == ProbabilityDistribution{Multivariate, PointMass}(Dict(:m=>[0.2]))
    @test ProbabilityDistribution(MatrixVariate, PointMass, m=transpose([0.2])) == ProbabilityDistribution{MatrixVariate, PointMass}(Dict(:m=>transpose([0.2])))
    @test ProbabilityDistribution(PointMass, m=0.2) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>0.2))
    @test ProbabilityDistribution(Univariate, PointMass) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0))
    @test ProbabilityDistribution(PointMass) == ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0))
    @test ProbabilityDistribution(Multivariate, PointMass) == ProbabilityDistribution{Multivariate, PointMass}(Dict(:m=>[1.0]))
    @test ProbabilityDistribution(MatrixVariate, PointMass) == ProbabilityDistribution{MatrixVariate, PointMass}(Dict(:m=>mat(1.0)))
    @test Message(PointMass) == Message{PointMass, Univariate}(ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0)))
    @test Message(Univariate, PointMass) == Message{PointMass, Univariate}(ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>1.0)))
    @test_throws Exception Message(Multivariate, PointMass, m=0.2)
end

@testset "ProbabilityDistribution and Message construction with functions" begin
    f(x) = x
    @test isa(ProbabilityDistribution(Univariate, Function, f=f).params[:f], Function)
    @test isa(ProbabilityDistribution(Univariate, Function, f=f), ProbabilityDistribution{Univariate, Function})
    @test isa(ProbabilityDistribution(Multivariate, Function, f=f), ProbabilityDistribution{Multivariate, Function})
    @test isa(ProbabilityDistribution(MatrixVariate, Function, f=f), ProbabilityDistribution{MatrixVariate, Function})
    @test isa(ProbabilityDistribution(Function, f=f), ProbabilityDistribution{Univariate, Function})
    @test isempty(vague(Function).params)
    @test isa(Message(Univariate, Function, f=()->()).dist.params[:f], Function)
    @test isa(Message(Univariate, Function), Message{Function, Univariate})
    @test isa(Message(Function), Message{Function, Univariate})
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

@testset "gaussianQuadrature" begin
    @test gaussianQuadrature(x -> 1.0, m=-1.0, v=2.0) == 0.9999996733487053
    @test gaussianQuadrature(x -> x, m=-1.0, v=2.0) == -0.9999996733487055
    @test gaussianQuadrature(x -> x^2, m=-1.0, v=2.0) == 3.0000013419448215
end

@testset "@RV" begin
    g = FactorGraph()

    # @RV should construct a new variable
    x = constant(0.0)
    @RV y ~ GaussianMeanVariance(x, constant(1.0))
    @test length(g.variables) == 3 # including constants
    @test y.id == :y # automatically assign id based on variable name in code
    @test haskey(g.variables, y.id)
    @test g.variables[y.id] === y

    # @RV ~ should reuse existing variable and handle keyword agruments
    y_old = y
    @RV y ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:g_node)
    @test length(g.variables) == 5 # including constants
    @test haskey(g.nodes, :g_node)
    @test y === y_old

    # @RV should handle array element assignments and explicit Variable ids
    g = FactorGraph()
    vars = Vector{Variable}(undef, 2)
    i = 1
    @RV [id=:v*i] vars[1] ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:tst1) # new Variable
    @test length(g.variables) == 3 # including constants
    @test vars[1].id == :v1
    @RV vars[1] ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:tst2) # existing Variable
    @test length(g.variables) == 5 # including constants
    @RV vars[2*i] ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @test vars[2*i].id == :vars_2
    varmatrix = Matrix{Variable}(undef,2,2)
    @RV varmatrix[1,2*i] ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @test varmatrix[1,2*i].id == :varmatrix_1_2
    vardict = Dict{Int,Variable}()
    @RV vardict[3*i+1] ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @test vardict[3*i+1].id == :vardict_4

    # @RV should work with '= syntax'
    g = FactorGraph()
    @RV x = constant(1.0) + constant(2.0)
    @test length(g.variables) == 3
    @test x.id == :x
    @RV [id=:my_y] y = x + constant(2.0)
    @test length(g.variables) == 5
    @test y.id == :my_y

    # @RV without node definition should create a new Variable
    g = FactorGraph()
    @RV x
    @test isa(x, Variable)
    @test g.variables[:x] === x
    @RV [id=:x_new] x
    @test g.variables[:x_new] === x
    @test length(g.variables) == 2
end

end #module