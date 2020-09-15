module NonlinearExtendedTest

using Test
using LinearAlgebra
using ForneyLab
using ForneyLab: outboundType, isApplicable, Extended, requiresBreaker, breakerParameters
using ForneyLab: SPNonlinearEOutNG, SPNonlinearEIn1GG, SPNonlinearEOutNGX, SPNonlinearEInGX, MNonlinearEInGX
using ForneyLab: concatenate, localLinearization, requiresBreaker, breakerParameters

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
    @test concatenate([[1.0, 2.0], [3.0]]) == ([1.0, 2.0, 3.0], [2, 1])
end

@testset "localLinearization" begin
    (a, b) = localLinearization(Univariate, g, 1.0)
    @test a == 2.0
    @test b == -6.0

    (A, b) = localLinearization(Multivariate, g, [1.0])
    @test A == mat(2.0)
    @test b == [-6.0]

    (A, b) = localLinearization(Univariate, h, [1.0, 2.0])
    @test A == [2.0, -1.0]'
    @test b == -1.0

    (A, b) = localLinearization(Multivariate, h, [[1.0], [2.0]])
    @test A == [2.0 -1.0]
    @test b == [-1.0]
end

@testset "requiresBreaker and breakerParameters" begin
    # Without given inverse
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Nonlinear{Extended}(y, x, g=g)
    
    @test requiresBreaker(nd.i[:out])
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())

    # With given inverse
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Nonlinear{Extended}(y, x, g=g, g_inv=g)
    
    @test !requiresBreaker(nd.i[:out])
    @test_throws Exception breakerParameters(nd.i[:out])
end


#-------------
# Update rules
#-------------

@testset "SPNonlinearEOutNG" begin
    @test SPNonlinearEOutNG <: SumProductRule{Nonlinear{Extended}}
    @test outboundType(SPNonlinearEOutNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearEOutNG, [Nothing, Message{GaussianMeanVariance}]) 

    @test ruleSPNonlinearEOutNG(g, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0)) == Message(Univariate, GaussianMeanVariance, m=-1.0, v=48.0)
    @test ruleSPNonlinearEOutNG(g, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[-1.0], v=mat(48.0))
end

@testset "SPNonlinearEOutNGX" begin
    @test SPNonlinearEOutNGX <: SumProductRule{Nonlinear{Extended}}
    @test outboundType(SPNonlinearEOutNGX) == Message{GaussianMeanVariance}
    @test !isApplicable(SPNonlinearEOutNGX, [Nothing, Message{Gaussian}]) 
    @test isApplicable(SPNonlinearEOutNGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test !isApplicable(SPNonlinearEOutNGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    @test ruleSPNonlinearEOutNGX(h, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=-1.0, v=49.0)
    @test ruleSPNonlinearEOutNGX(h, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[-1.0], v=mat(49.0))
end

@testset "SPNonlinearEIn1GG" begin
    @test SPNonlinearEIn1GG <: SumProductRule{Nonlinear{Extended}}
    @test outboundType(SPNonlinearEIn1GG) == Message{Gaussian}
    @test isApplicable(SPNonlinearEIn1GG, [Message{Gaussian}, Nothing]) 

    # Without given inverse
    @test ruleSPNonlinearEIn1GG(g, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=14.666666666666666, w=5.333333333333333)
    @test ruleSPNonlinearEIn1GG(g, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[14.666666666666671], w=mat(5.333333333333335))

    # With given inverse
    @test ruleSPNonlinearEIn1GG(g, g_inv, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.6457513110645907, v=0.10714285714285711)
    @test ruleSPNonlinearEIn1GG(g, g_inv, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.6457513110645907], v=mat(0.10714285714285711))
end

@testset "SPNonlinearEInGX" begin
    @test SPNonlinearEInGX <: SumProductRule{Nonlinear{Extended}}
    @test outboundType(SPNonlinearEInGX) == Message{Gaussian}
    @test !isApplicable(SPNonlinearEInGX, [Message{Gaussian}, Nothing]) 
    @test !isApplicable(SPNonlinearEInGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test isApplicable(SPNonlinearEInGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    # Without given inverse
    @test ruleSPNonlinearEInGX(h, 1, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=10.999999999999996, w=3.9999999999999982)
    @test ruleSPNonlinearEInGX(h, 1, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[10.999999999999996], w=mat(3.9999999999999982))

    # With given inverse
    @test ruleSPNonlinearEInGX(h, h_inv_x, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=1.3228756555322954, v=0.14285714285714282)
    @test ruleSPNonlinearEInGX(h, h_inv_x, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.3228756555322954], v=mat(0.14285714285714282))
end

@testset "MNonlinearEInGX" begin
    @test MNonlinearEInGX <: MarginalRule{Nonlinear{Extended}}
    @test isApplicable(MNonlinearEInGX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleMNonlinearEInGX(h, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.6, 4.85], v=[0.20000000000000007 0.19999999999999998; 0.19999999999999998 0.95])
    @test ruleMNonlinearEInGX(h, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.6, 4.85], v=[0.20000000000000007 0.19999999999999998; 0.19999999999999998 0.95])
end


#------------
# Integration
#------------

@testset "Nonlinear integration via local linear approximation with given inverse" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{Extended}(y, x, g=g, g_inv=g_inv)
    
    @test isa(n, Nonlinear{Extended})
    
    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEOutNG(g, nothing, messages[2])", code)
    @test !occursin("g_inv", code)

    # Backward; g_inv should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEIn1GG(g, g_inv, messages[2], nothing)", code)
end

@testset "Multi-argument nonlinear integration via local linear approximation" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    @RV z ~ GaussianMeanVariance(5.0, 1.0)
    n = Nonlinear{Extended}(y, x, z, g=h, g_inv=[h_inv_x, nothing])
    
    # Forward; h_inv_x should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEOutNGX(h, nothing, messages[3], messages[1])", code)
    @test !occursin("h_inv_x", code)

    # Backward with given inverse; h_inv_x should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEInGX(h, h_inv_x, messages[3], nothing, messages[1])", code)

    # Backward without given inverse
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(z)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEInGX(h, 2, messages[3], messages[2], messages[1])", code)
    @test !occursin("h_inv_x", code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", code)
end

@testset "Nonlinear integration via local linear approximation without given inverse" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{Extended}(y, x, g=g)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEOutNG(g, nothing, messages[1])", code)
    @test !occursin("g_inv", code)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearEIn1GG(g, messages[2], messages[1])", code)
    @test !occursin("g_inv", code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", code)
end

end #module