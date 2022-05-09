module DeltaExtendedTest

using Test
using LinearAlgebra
using ForneyLab
using ForneyLab: outboundType, isApplicable, Extended, requiresBreaker, breakerParameters
using ForneyLab: SPDeltaEOutNG, SPDeltaEIn1GG, SPDeltaEOutNGX, SPDeltaEInGX, MDeltaEInGX
using ForneyLab: concatenate, localLinearizationSingleIn, localLinearizationMultiIn, requiresBreaker, breakerParameters

f(x) = x

g(x::Number) = x^2 - 5.0
g(x::Vector) = x.^2 .- 5.0
g_inv(y::Number) = sqrt(y + 5.0)
g_inv(y::Vector) = sqrt.(y .+ 5.0)

h(x::Number, y::Number) = x^2 - y
h(x::Vector, y::Vector) = x.^2 .- y
h_inv_x(z::Number, y::Number) = sqrt(z + y)
h_inv_x(z::Vector, y::Vector) = sqrt.(z .+ y)


#--------
# Helpers
#--------

@testset "concatenate" begin
    @test concatenate([[1.0, 2.0], [3.0]]) == ([1.0, 2.0, 3.0], [(2,), (1,)])
    @test concatenate([[1.0, 2.0], 3.0]) == ([1.0, 2.0, 3.0], [(2,), ()])
    @test concatenate([[1.0, 2.0], [3.0 5.0; 4.0 6.0]]) == ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [(2,), (2,2)])
end

@testset "localLinearization" begin
    (a, b) = localLinearizationSingleIn(g, 1.0)
    @test a == 2.0
    @test b == -6.0

    (A, b) = localLinearizationSingleIn(g, [1.0])
    @test A == mat(2.0)
    @test b == [-6.0]

    (A, b) = localLinearizationMultiIn(h, [1.0, 2.0])
    @test A == [2.0, -1.0]'
    @test b == -1.0

    (A, b) = localLinearizationMultiIn(h, [[1.0], [2.0]])
    @test A == [2.0 -1.0]
    @test b == [-1.0]
end

@testset "Delta construction alias" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = Nonlinear{Extended}(y, x, g=g)
    
    @test typeof(nd) == Delta{Extended}
    @test nd.id == :delta_1
end

@testset "requiresBreaker and breakerParameters" begin
    # Without given inverse
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = Gaussian{Moments}(x, 0.0, 1.0)
    Delta{Extended}(y, x, g=g)
    
    @test requiresBreaker(nd.i[:out])
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{Gaussian{Moments}, Univariate}, ())

    # With given inverse
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = Gaussian{Moments}(x, 0.0, 1.0)
    Delta{Extended}(y, x, g=g, g_inv=g)
    
    @test !requiresBreaker(nd.i[:out])
    @test_throws Exception breakerParameters(nd.i[:out])
end


#-------------
# Update rules
#-------------

@testset "SPDeltaEOutNG" begin
    @test SPDeltaEOutNG <: SumProductRule{Delta{Extended}}
    @test outboundType(SPDeltaEOutNG) == Message{Gaussian{Moments}}
    @test isApplicable(SPDeltaEOutNG, [Nothing, Message{Gaussian{Moments}}]) 

    @test ruleSPDeltaEOutNG(g, nothing, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0)) == Message(Univariate, Gaussian{Moments}, m=-1.0, v=48.0)
    @test ruleSPDeltaEOutNG(g, nothing, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0))) == Message(Multivariate, Gaussian{Moments}, m=[-1.0], v=mat(48.0))
end

@testset "SPDeltaEOutNGX" begin
    @test SPDeltaEOutNGX <: SumProductRule{Delta{Extended}}
    @test outboundType(SPDeltaEOutNGX) == Message{Gaussian{Moments}}
    @test !isApplicable(SPDeltaEOutNGX, [Nothing, Message{Gaussian}]) 
    @test isApplicable(SPDeltaEOutNGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test !isApplicable(SPDeltaEOutNGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    @test ruleSPDeltaEOutNGX(h, nothing, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0), Message(Univariate, Gaussian{Moments}, m=5.0, v=1.0)) == Message(Univariate, Gaussian{Moments}, m=-1.0, v=49.0)
    @test ruleSPDeltaEOutNGX(h, nothing, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0)), Message(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(1.0))) == Message(Multivariate, Gaussian{Moments}, m=[-1.0], v=mat(49.0))
end

@testset "SPDeltaEIn1GG" begin
    @test SPDeltaEIn1GG <: SumProductRule{Delta{Extended}}
    @test outboundType(SPDeltaEIn1GG) == Message{Gaussian}
    @test isApplicable(SPDeltaEIn1GG, [Message{Gaussian}, Nothing]) 

    # Without given inverse
    @test ruleSPDeltaEIn1GG(g, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0), Message(Univariate, Gaussian{Moments}, m=2.0, v=1.0)) == Message(Univariate, Gaussian{Canonical}, xi=14.666666666666666, w=5.333333333333333)
    @test ruleSPDeltaEIn1GG(g, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0)), Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(1.0))) == Message(Multivariate, Gaussian{Canonical}, xi=[14.666666666666671], w=mat(5.333333333333335))

    # With given inverse
    @test ruleSPDeltaEIn1GG(g, g_inv, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0), nothing) == Message(Univariate, Gaussian{Moments}, m=2.6457513110645907, v=0.10714285714285711)
    @test ruleSPDeltaEIn1GG(g, g_inv, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0)), nothing) == Message(Multivariate, Gaussian{Moments}, m=[2.6457513110645907], v=mat(0.10714285714285711))
end

@testset "SPDeltaEInGX" begin
    @test SPDeltaEInGX <: SumProductRule{Delta{Extended}}
    @test outboundType(SPDeltaEInGX) == Message{Gaussian}
    @test !isApplicable(SPDeltaEInGX, [Message{Gaussian}, Nothing]) 
    @test !isApplicable(SPDeltaEInGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test isApplicable(SPDeltaEInGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    # Without given inverse
    @test ruleSPDeltaEInGX(h, 1, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0), Message(Univariate, Gaussian{Moments}, m=2.0, v=1.0), Message(Univariate, Gaussian{Moments}, m=5.0, v=1.0)) == Message(Univariate, Gaussian{Canonical}, xi=10.999999999999996, w=3.9999999999999982)
    @test ruleSPDeltaEInGX(h, 1, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0)), Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(1.0)), Message(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(1.0))) == Message(Multivariate, Gaussian{Canonical}, xi=[10.999999999999996], w=mat(3.9999999999999982))

    # With given inverse
    @test ruleSPDeltaEInGX(h, h_inv_x, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0), nothing, Message(Univariate, Gaussian{Moments}, m=5.0, v=1.0)) == Message(Univariate, Gaussian{Moments}, m=2.6457513110645907, v=0.14285714285714282)
    @test ruleSPDeltaEInGX(h, h_inv_x, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0)), nothing, Message(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(1.0))) == Message(Multivariate, Gaussian{Moments}, m=[2.6457513110645907], v=mat(0.14285714285714282))
end

@testset "MDeltaEInGX" begin
    @test MDeltaEInGX <: MarginalRule{Delta{Extended}}
    @test isApplicable(MDeltaEInGX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleMDeltaEInGX(h, Message(Univariate, Gaussian{Moments}, m=2.0, v=3.0), Message(Univariate, Gaussian{Moments}, m=2.0, v=1.0), Message(Univariate, Gaussian{Moments}, m=5.0, v=1.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[2.6, 4.85], v=[0.20000000000000007 0.19999999999999998; 0.19999999999999998 0.95])
    @test ruleMDeltaEInGX(h, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(3.0)), Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(1.0)), Message(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(1.0))) == Distribution(Multivariate, Gaussian{Moments}, m=[2.6, 4.85], v=[0.20000000000000007 0.19999999999999998; 0.19999999999999998 0.95])
end


#------------
# Integration
#------------

@testset "Delta integration via local linear approximation with given inverse" begin
    fg = FactorGraph()

    @RV x ~ Gaussian{Moments}(2.0, 1.0)
    @RV y ~ Gaussian{Moments}(2.0, 3.0)
    n = Delta{Extended}(y, x, g=g, g_inv=g_inv)
    
    @test isa(n, Delta{Extended})
    
    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEOutNG(g, nothing, messages[2])", code)
    @test !occursin("g_inv", code)

    # Backward; g_inv should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEIn1GG(g, g_inv, messages[2], nothing)", code)
end

@testset "Multi-argument integration via local linear approximation" begin
    fg = FactorGraph()

    @RV x ~ Gaussian{Moments}(2.0, 1.0)
    @RV y ~ Gaussian{Moments}(2.0, 3.0)
    @RV z ~ Gaussian{Moments}(5.0, 1.0)
    n = Delta{Extended}(y, x, z, g=h, g_inv=[h_inv_x, nothing])
    
    # Forward; h_inv_x should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEOutNGX(h, nothing, messages[3], messages[1])", code)
    @test !occursin("h_inv_x", code)

    # Backward with given inverse; h_inv_x should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEInGX(h, h_inv_x, messages[3], nothing, messages[1])", code)

    # Backward without given inverse
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(z)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEInGX(h, 2, messages[3], messages[2], messages[1])", code)
    @test !occursin("h_inv_x", code)
    @test occursin("messages[1] = Message(vague(Gaussian{Moments}))", code)
end

@testset "Delta integration via local linear approximation without given inverse" begin
    fg = FactorGraph()

    @RV x ~ Gaussian{Moments}(2.0, 1.0)
    @RV y ~ Gaussian{Moments}(2.0, 3.0)
    n = Delta{Extended}(y, x, g=g)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEOutNG(g, nothing, messages[1])", code)
    @test !occursin("g_inv", code)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaEIn1GG(g, messages[2], messages[1])", code)
    @test !occursin("g_inv", code)
    @test occursin("messages[1] = Message(vague(Gaussian{Moments}))", code)
end

end #module