module DeltaSamplingTest

using Test
using Random
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, logPdf, unsafeMean, unsafeVar, Sampling, requiresBreaker, breakerParameters
using ForneyLab: SPDeltaSInGX, SPDeltaSOutNM, SPDeltaSIn1MN, SPDeltaSOutNMX, SPDeltaSInMX, MDeltaSInMGX, gradientOptimization
using ForwardDiff

Random.seed!(1234)

f(x) = x
g(x) = x^2
h(x, y) = x + y

@testset "gradientOptimization" begin
    f(x) = -x[1]^2 - x[2]^2
    grad(x) = ForwardDiff.gradient(f, x)
    res = gradientOptimization(f, grad, [4.0, 1.0], 0.001)

    # Optimum is [0.0, 0.0]
    @test isapprox(res, [0.0, 0.0], atol=1e-16)
end

@testset "requiresBreaker and breakerParameters" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Delta{Sampling}(y, x, g=g)

    @test !requiresBreaker(nd.i[:out]) # Single-input Delta{Sampling} does not require breaker

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    y = Variable()
    GaussianMeanVariance(z, 0.0, 1.0)
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Delta{Sampling}(y, x, z, g=h)

    @test requiresBreaker(nd.i[:out]) # Multi-input Delta{Sampling} requires breaker
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())
end


#-------------
# Update rules
#-------------

@testset "SPDeltaSOutNM" begin
    @test SPDeltaSOutNM <: SumProductRule{Delta{Sampling}}
    @test outboundType(SPDeltaSOutNM) == Message{SampleList}
    @test isApplicable(SPDeltaSOutNM, [Nothing, Message{Gaussian}])
    @test isApplicable(SPDeltaSOutNM, [Nothing, Message{Bernoulli}])

    msg = ruleSPDeltaSOutNM(f, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny), n_samples=1)
    @test isapprox(msg.dist.params[:s][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
    msg = ruleSPDeltaSOutNM(f, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(tiny)), n_samples=1)
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "SPDeltaSIn1MN" begin
    @test SPDeltaSIn1MN <: SumProductRule{Delta{Sampling}}
    @test outboundType(SPDeltaSIn1MN) == Message{Function}
    @test isApplicable(SPDeltaSIn1MN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])

    log_pdf(x) = ruleSPDeltaSIn1MN(f, Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), nothing, n_samples=0).dist.params[:log_pdf](x)
    @test log_pdf(1.5) == logPdf(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), 1.5)
end

@testset "SPDeltaSInGX" begin
    @test SPDeltaSInGX <: SumProductRule{Delta{Sampling}}
    @test outboundType(SPDeltaSInGX) == Message{GaussianWeightedMeanPrecision}
    @test !isApplicable(SPDeltaSInGX, [Nothing, Message{Gamma}])
    @test isApplicable(SPDeltaSInGX, [Message{Gaussian}, Message{Gaussian}, Nothing])

    msg_out = Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    msg_in2 = Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0)

    res = ruleSPDeltaSInGX(h, 1, msg_out, msg_in1, msg_in2, n_samples=1000)

    @test isapprox(mean(res.dist), 0.9998084187686186, atol=0.1)
    @test isapprox(var(res.dist), 1.9999999999999982, atol=0.1)
end

@testset "prod!" begin
    d = prod!(Distribution(Multivariate, Function, log_pdf=(s)->s), Distribution(Multivariate, Function, log_pdf=(s)->s))
    @test isa(d, Distribution{Multivariate,Function})
    @test d.params[:log_pdf](1) == 2
end

@testset "SPDeltaSOutNMX" begin
    @test SPDeltaSOutNMX <: SumProductRule{Delta{Sampling}}
    @test outboundType(SPDeltaSOutNMX) == Message{SampleList}
    @test !isApplicable(SPDeltaSOutNMX, [Nothing, Message{Gaussian}])
    @test isApplicable(SPDeltaSOutNMX, [Nothing, Message{Gaussian}, Message{Gamma}])
    @test isApplicable(SPDeltaSOutNMX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    msg = ruleSPDeltaSOutNMX(h, nothing, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0), n_samples=1)
    @test typeof(msg) == Message{SampleList, Univariate}
end

@testset "SPDeltaSInMX" begin
    @test SPDeltaSInMX <: SumProductRule{Delta{Sampling}}
    @test outboundType(SPDeltaSInMX) == Message{Function}
    @test !isApplicable(SPDeltaSInMX, [Message, Message{Gaussian}])
    @test isApplicable(SPDeltaSInMX, [Message, Nothing, Message{Gaussian}, Message{Gamma}])
    @test isApplicable(SPDeltaSInMX, [Message, Nothing, Message{PointMass}])
    @test isApplicable(SPDeltaSInMX, [Message, Message{PointMass}, Nothing])

    msg = ruleSPDeltaSInMX(h, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), nothing, Message(Univariate, PointMass, m=2.0))
    @test msg.dist.params[:log_pdf](1.0) == 0.23235401329235006
    msg = ruleSPDeltaSInMX(h, 1, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0))
    @test typeof(msg.dist.params[:log_pdf](1.0)) == Float64
end

@testset "MDeltaSInMGX" begin
    @test MDeltaSInMGX <: MarginalRule{Delta{Sampling}}
    @test isApplicable(MDeltaSInMGX, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(MDeltaSInMGX, [Nothing, Message{Gaussian}])
    @test !isApplicable(MDeltaSInMGX, [Nothing, Message{Gaussian}, Message{Gamma}])

    dist = ruleMDeltaSInMGX(h, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0))
    @test dist == Distribution(Multivariate, GaussianMeanPrecision, m=[1.0, 2.0], w=[11.0 10.0; 10.0 10.5])
end


#------------
# Integration
#------------

@testset "Delta integration via sampling" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Delta{Sampling}(y, x, g=g, n_samples=2000)

    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaSOutNM(g, nothing, messages[2], n_samples=2000)", code)
end

@testset "Delta integration via sampling" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    @RV z ~ GaussianMeanVariance(2.0, 3.0)
    n = Delta{Sampling}(z, x, y, g=g)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)

    @test occursin("ruleSPDeltaSInGX(g, 2, messages[3], messages[2], messages[1])", code)
end

@testset "Delta integration via sampling with specified variate types" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance([2.0], mat(1.0))
    @RV y ~ GaussianMeanVariance([2.0], mat(3.0))
    @RV z ~ GaussianMeanVariance([2.0], mat(3.0))
    n = Delta{Sampling}(z, x, y, g=g, dims=[(1,), (1,), (1,)])

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    
    @test occursin("ruleSPDeltaSInGX(g, 2, messages[3], messages[2], messages[1], dims=(1,))", code)
end

end # module