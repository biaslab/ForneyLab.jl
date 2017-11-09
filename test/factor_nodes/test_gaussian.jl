module GaussianWeightedMeanPrecisionTest

using Base.Test
using ForneyLab
import ForneyLab: ensureParameters!, prod!, ==, unsafeMean, unsafeVar, unsafeCov, isWellDefined, isProper, sample, dims

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, Gaussian, m=ones(2), v=diageye(2))) == 2
end

@testset "vague" begin
    # Univariate
    @test vague(Gaussian) == ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=huge)

    # Multivariate
    @test vague(Gaussian, 2) == ProbabilityDistribution(Multivariate, Gaussian, m=[0.0, 0.0], v=huge*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0))
    @test !isProper(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=-1.0))
    @test !isProper(ProbabilityDistribution(Univariate, Gaussian, m=0.0, w=-1.0))

    # Multivariate
    @test isProper(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].'))
    @test isProper(ProbabilityDistribution(Multivariate, Gaussian, m=ones(2), v=diageye(2)))
    @test !isProper(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[-1.0].'))
    @test !isProper(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], w=[-1.0].'))
end

@testset "sample" begin
    # Univariate
    @test isa(sample(ProbabilityDistribution(Univariate, Gaussian, m=1.0, v=2.0)), Float64)

    # Multivariate
    @test isa(sample(ProbabilityDistribution(Multivariate, Gaussian, m=[1.2, 2.7], v=[2.0 -0.5; -0.5 1.5])), Vector{Float64})
    @test isa(sample(ProbabilityDistribution(Multivariate, Gaussian, m=[1.2, 2.7], v=Diagonal([2.0, 1.5]))), Vector{Float64})
end

@testset "prod!" begin
    # Univariate
    @test ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0) * ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0) == ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=2.0)
    @test ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0) * ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0) == ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=2.0)

    # Multivariate
    @test ProbabilityDistribution(Multivariate, Gaussian, xi=zeros(2), w=diageye(2)) * ProbabilityDistribution(Multivariate, Gaussian, xi=zeros(2), w=diageye(2)) == ProbabilityDistribution(Multivariate, Gaussian, xi=zeros(2), w=2.0*diageye(2))
    @test ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), v=diageye(2)) * ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), v=diageye(2)) == ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), v=0.5*diageye(2))
end

@testset "isWellDefined" begin
    # Univariate
    @test isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0))
    @test isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=0.0, w=1.0))
    @test isWellDefined(ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0))
    @test isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=0.0, xi=0.0, w=1.0, v=1.0))
    @test !isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=NaN, v=1.0))
    @test !isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=0.0, w=NaN))
    @test !isWellDefined(ProbabilityDistribution(Univariate, Gaussian, v=1.0, w=1.0))
    @test !isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=0.0, xi=0.0))
    @test !isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=NaN, v=1.0, w=1.0))
    @test !isWellDefined(ProbabilityDistribution(Univariate, Gaussian, m=0.0, xi=0.0, v=NaN))

    # Multivariate
    @test isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].'))
    @test isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], w=[1.0].'))
    @test isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), w=diageye(2)))
    @test isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, xi=[0.0], w=[1.0].'))
    @test isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], xi=[0.0], w=[1.0].', v=[1.0].'))
    @test !isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[NaN], v=[1.0].'))
    @test !isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], w=[NaN].'))
    @test !isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, v=[1.0].', w=[1.0].'))
    @test !isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], xi=[0.0]))
    @test !isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[NaN], v=[1.0].', w=[1.0].'))
    @test !isWellDefined(ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], xi=[0.0], v=[NaN].'))
end

@testset "ensureParameters!" begin
    # Univariate
    dist = ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)
    @test ensureParameters!(dist, (:m, :v)).params == Dict(:m=>0.0, :v=>1.0)
    dist = ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)
    @test ensureParameters!(dist, (:xi, :w)).params == Dict(:m=>0.0, :v=>1.0, :xi=>0.0, :w=>1.0)

    # Multivariate
    dist = ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].')
    @test ensureParameters!(dist, (:m, :v)).params == Dict(:m=>[0.0], :v=>[1.0].')
    dist = ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].')
    @test ensureParameters!(dist, (:xi, :w)).params == Dict(:m=>[0.0], :v=>[1.0].', :xi=>[0.0], :w=>[1.0].')
end

@testset "==" begin
    # Univariate
    @test ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0) == ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0)
    @test ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0) == ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)

    # Multivariate
    @test ProbabilityDistribution(Multivariate, Gaussian, xi=[0.0], w=[1.0].') == ProbabilityDistribution(Multivariate, Gaussian, xi=[0.0], w=[1.0].')
    @test ProbabilityDistribution(Multivariate, Gaussian, xi=[0.0], w=[1.0].') == ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].')
end

@testset "unsafe mean and variance" begin
    # Univariate
    @test unsafeMean(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=2.0)) == 0.0
    @test unsafeVar(ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=2.0)) == 2.0
    @test unsafeCov(ProbabilityDistribution(Univariate, Gaussian, m=0.0, w=2.0)) == 0.5

    # Multivariate
    @test unsafeMean(ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), v=2.0*eye(2))) == zeros(2)
    @test unsafeVar(ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), v=2.0*eye(2))) == 2.0*ones(2)
    @test unsafeCov(ProbabilityDistribution(Multivariate, Gaussian, m=zeros(2), w=2.0*eye(2))) == [0.4999999999999999 -0.0; -0.0 0.4999999999999999]
end

end #module