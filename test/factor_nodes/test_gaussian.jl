module GaussianWeightedMeanPrecisionTest

using Base.Test
using ForneyLab
import ForneyLab: ensureParameters!, prod!, ==, unsafeMean, unsafeVar

@testset "vague" begin
    # Univariate
    @test vague(Univariate{Gaussian}) == Univariate(Gaussian, m=0.0, v=huge)

    # Multivariate
    @test vague(Multivariate{Gaussian, 2}) == Multivariate(Gaussian, m=[0.0, 0.0], V=huge*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Univariate(Gaussian, m=0.0, v=1.0))
    @test !isProper(Univariate(Gaussian, m=0.0, v=-1.0))
    @test !isProper(Univariate(Gaussian, m=0.0, w=-1.0))

    # Multivariate
    @test isProper(Multivariate(Gaussian, m=[0.0], V=[1.0].'))
    @test isProper(Multivariate(Gaussian, m=ones(2), V=diageye(2)))
    @test !isProper(Multivariate(Gaussian, m=[0.0], V=[-1.0].'))
    @test !isProper(Multivariate(Gaussian, m=[0.0], W=[-1.0].'))
end

@testset "sample" begin
    # Univariate
    @test isa(Float64, sample(Univariate(Gaussian, m=1.0, v=2.0)))

    # Multivariate
    @test isa(Vector{Float64}, sample(Multivariate(Gaussian, m=[1.2, 2.7], V=[2.0 -0.5; -0.5 1.5])))
    @test isa(Vector{Float64}, sample(Multivariate(Gaussian, m=[1.2, 2.7], V=Diagonal([2.0, 1.5]))))
end

@testset "prod!" begin
    # Univariate
    @test Univariate(Gaussian, xi=0.0, w=1.0) * Univariate(Gaussian, xi=0.0, w=1.0) == Univariate(Gaussian, xi=0.0, w=2.0)
    @test Univariate(Gaussian, m=0.0, v=1.0) * Univariate(Gaussian, m=0.0, v=1.0) == Univariate(Gaussian, xi=0.0, w=2.0)

    # Multivariate
    @test Multivariate(Gaussian, xi=zeros(2), W=diageye(2)) * Multivariate(Gaussian, xi=zeros(2), W=diageye(2)) == Multivariate(Gaussian, xi=zeros(2), W=2.0*diageye(2))
    @test Multivariate(Gaussian, m=zeros(2), V=diageye(2)) * Multivariate(Gaussian, m=zeros(2), V=diageye(2)) == Multivariate(Gaussian, m=zeros(2), V=0.5*diageye(2))
end

@testset "isWellDefined" begin
    # Univariate
    @test isWellDefined(Univariate(Gaussian, m=0.0, v=1.0))
    @test isWellDefined(Univariate(Gaussian, m=0.0, w=1.0))
    @test isWellDefined(Univariate(Gaussian, xi=0.0, w=1.0))
    @test isWellDefined(Univariate(Gaussian, m=0.0, xi=0.0, w=1.0, v=1.0))
    @test !isWellDefined(Univariate(Gaussian, m=NaN, v=1.0))
    @test !isWellDefined(Univariate(Gaussian, m=0.0, w=NaN))
    @test !isWellDefined(Univariate(Gaussian, v=1.0, w=1.0))
    @test !isWellDefined(Univariate(Gaussian, m=0.0, xi=0.0))
    @test !isWellDefined(Univariate(Gaussian, m=NaN, v=1.0, w=1.0))
    @test !isWellDefined(Univariate(Gaussian, m=0.0, xi=0.0, v=NaN))

    # Multivariate
    @test isWellDefined(Multivariate(Gaussian, m=[0.0], V=[1.0].'))
    @test isWellDefined(Multivariate(Gaussian, m=[0.0], W=[1.0].'))
    @test isWellDefined(Multivariate(Gaussian, m=zeros(2), W=diageye(2)))
    @test isWellDefined(Multivariate(Gaussian, xi=[0.0], W=[1.0].'))
    @test isWellDefined(Multivariate(Gaussian, m=[0.0], xi=[0.0], W=[1.0].', V=[1.0].'))
    @test !isWellDefined(Multivariate(Gaussian, m=[NaN], V=[1.0].'))
    @test !isWellDefined(Multivariate(Gaussian, m=[0.0], W=[NaN].'))
    @test !isWellDefined(Multivariate(Gaussian, V=[1.0].', W=[1.0].'))
    @test !isWellDefined(Multivariate(Gaussian, m=[0.0], xi=[0.0]))
    @test !isWellDefined(Multivariate(Gaussian, m=[NaN], V=[1.0].', W=[1.0].'))
    @test !isWellDefined(Multivariate(Gaussian, m=[0.0], xi=[0.0], V=[NaN].'))
end

@testset "ensureParameters!" begin
    # Univariate
    dist = Univariate(Gaussian, m=0.0, v=1.0)
    @test ensureParameters!(dist, (:m, :v)).params == Dict(:m=>0.0, :v=>1.0)
    dist = Univariate(Gaussian, m=0.0, v=1.0)
    @test ensureParameters!(dist, (:xi, :w)).params == Dict(:m=>0.0, :v=>1.0, :xi=>0.0, :w=>1.0)

    # Multivariate
    dist = Multivariate(Gaussian, m=[0.0], V=[1.0].')
    @test ensureParameters!(dist, (:m, :V)).params == Dict(:m=>[0.0], :V=>[1.0].')
    dist = Multivariate(Gaussian, m=[0.0], V=[1.0].')
    @test ensureParameters!(dist, (:xi, :W)).params == Dict(:m=>[0.0], :V=>[1.0].', :xi=>[0.0], :W=>[1.0].')
end

@testset "==" begin
    # Univariate
    @test Univariate(Gaussian, xi=0.0, w=1.0) == Univariate(Gaussian, xi=0.0, w=1.0)
    @test Univariate(Gaussian, xi=0.0, w=1.0) == Univariate(Gaussian, m=0.0, v=1.0)

    # Multivariate
    @test Multivariate(Gaussian, xi=[0.0], W=[1.0].') == Multivariate(Gaussian, xi=[0.0], W=[1.0].')
    @test Multivariate(Gaussian, xi=[0.0], W=[1.0].') == Multivariate(Gaussian, m=[0.0], V=[1.0].')
end

@testset "unsafe mean and variance" begin
    # Univariate
    @test unsafeMean(Univariate(Gaussian, m=0.0, v=2.0)) == 0.0
    @test unsafeVar(Univariate(Gaussian, m=0.0, v=2.0)) == 2.0
    @test unsafeCov(Univariate(Gaussian, m=0.0, w=2.0)) == 0.5

    # Multivariate
    @test unsafeMean(Multivariate(Gaussian, m=zeros(2), v=2.0*eye(2))) == zeros(2)
    @test unsafeVar(Multivariate(Gaussian, m=zeros(2), v=2.0*eye(2))) == 2.0*ones(2)
    @test unsafeCov(Multivariate(Gaussian, m=zeros(2), w=2.0*eye(2))) == 0.5*eye(2)
end





end #module