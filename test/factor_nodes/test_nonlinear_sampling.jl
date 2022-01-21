module NonlinearSamplingTest

using Test
using Random
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, logPdf, unsafeMean, unsafeVar, Sampling, requiresBreaker, breakerParameters
using ForneyLab: SPNonlinearSInGX, SPNonlinearSOutNM, SPNonlinearSIn1MN, SPNonlinearSOutNMX, SPNonlinearSInMX, MNonlinearSInMGX, gradientOptimization
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
    Nonlinear{Sampling}(y, x, g=g)

    @test !requiresBreaker(nd.i[:out]) # Single-input Nonlinear{Sampling} does not require breaker

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    y = Variable()
    GaussianMeanVariance(z, 0.0, 1.0)
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Nonlinear{Sampling}(y, x, z, g=h)

    @test requiresBreaker(nd.i[:out]) # Multi-input Nonlinear{Sampling} requires breaker
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())
end


#-------------
# Update rules
#-------------

@testset "SPNonlinearSOutNM" begin
    @test SPNonlinearSOutNM <: SumProductRule{Nonlinear{Sampling}}
    @test outboundType(SPNonlinearSOutNM) == Message{SampleList}
    @test isApplicable(SPNonlinearSOutNM, [Nothing, Message{Gaussian}])
    @test isApplicable(SPNonlinearSOutNM, [Nothing, Message{Bernoulli}])

    msg = ruleSPNonlinearSOutNM(f, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny), n_samples=1)
    @test isapprox(msg.dist.params[:s][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
    msg = ruleSPNonlinearSOutNM(f, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(tiny)), n_samples=1)
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "SPNonlinearSIn1MN" begin
    @test SPNonlinearSIn1MN <: SumProductRule{Nonlinear{Sampling}}
    @test outboundType(SPNonlinearSIn1MN) == Message{Function}
    @test isApplicable(SPNonlinearSIn1MN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])

    log_pdf(x) = ruleSPNonlinearSIn1MN(f, Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), nothing, n_samples=0).dist.params[:log_pdf](x)
    @test log_pdf(1.5) == logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), 1.5)
end

@testset "SPNonlinearSInGX" begin
    @test SPNonlinearSInGX <: SumProductRule{Nonlinear{Sampling}}
    @test outboundType(SPNonlinearSInGX) == Message{GaussianWeightedMeanPrecision}
    @test !isApplicable(SPNonlinearSInGX, [Nothing, Message{Gamma}])
    @test isApplicable(SPNonlinearSInGX, [Message{Gaussian}, Message{Gaussian}, Nothing])

    msg_out = Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    msg_in2 = Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0)

    res = ruleSPNonlinearSInGX(h, 1, msg_out, msg_in1, msg_in2, n_samples=1000)

    @test isapprox(mean(res.dist), 0.9998084187686186, atol=0.1)
    @test isapprox(var(res.dist), 1.9999999999999982, atol=0.1)
end

@testset "prod!" begin
    d = prod!(ProbabilityDistribution(Multivariate, Function, log_pdf=(s)->s), ProbabilityDistribution(Multivariate, Function, log_pdf=(s)->s))
    @test isa(d, ProbabilityDistribution{Multivariate,Function})
    @test d.params[:log_pdf](1) == 2
end

@testset "SPNonlinearSOutNMX" begin
    @test SPNonlinearSOutNMX <: SumProductRule{Nonlinear{Sampling}}
    @test outboundType(SPNonlinearSOutNMX) == Message{SampleList}
    @test !isApplicable(SPNonlinearSOutNMX, [Nothing, Message{Gaussian}])
    @test isApplicable(SPNonlinearSOutNMX, [Nothing, Message{Gaussian}, Message{Gamma}])
    @test isApplicable(SPNonlinearSOutNMX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    msg = ruleSPNonlinearSOutNMX(h, nothing, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0), n_samples=1)
    @test typeof(msg) == Message{SampleList, Univariate}
end

@testset "SPNonlinearSInMX" begin
    @test SPNonlinearSInMX <: SumProductRule{Nonlinear{Sampling}}
    @test outboundType(SPNonlinearSInMX) == Message{Function}
    @test !isApplicable(SPNonlinearSInMX, [Message, Message{Gaussian}])
    @test isApplicable(SPNonlinearSInMX, [Message, Nothing, Message{Gaussian}, Message{Gamma}])
    @test isApplicable(SPNonlinearSInMX, [Message, Nothing, Message{PointMass}])
    @test isApplicable(SPNonlinearSInMX, [Message, Message{PointMass}, Nothing])

    msg = ruleSPNonlinearSInMX(h, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), nothing, Message(Univariate, PointMass, m=2.0))
    @test msg.dist.params[:log_pdf](1.0) == 0.23235401329235006
    msg = ruleSPNonlinearSInMX(h, 1, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0))
    @test typeof(msg.dist.params[:log_pdf](1.0)) == Float64
end

@testset "MNonlinearSInMGX" begin
    @test MNonlinearSInMGX <: MarginalRule{Nonlinear{Sampling}}
    @test isApplicable(MNonlinearSInMGX, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(MNonlinearSInMGX, [Nothing, Message{Gaussian}])
    @test !isApplicable(MNonlinearSInMGX, [Nothing, Message{Gaussian}, Message{Gamma}])

    dist = ruleMNonlinearSInMGX(h, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0))
    @test dist == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[1.0, 2.0], w=[11.0 10.0; 10.0 10.5])
end


#------------
# Integration
#------------

@testset "Nonlinear integration via sampling" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{Sampling}(y, x, g=g, n_samples=2000)

    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearSOutNM(g, nothing, messages[2], n_samples=2000)", code)
end

@testset "Nonlinear integration via sampling" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    @RV z ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{Sampling}(z, x, y, g=g)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)

    @test occursin("ruleSPNonlinearSInGX(g, 2, messages[3], messages[2], messages[1])", code)
end

@testset "Nonlinear integration via sampling with specified variate types" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance([2.0], mat(1.0))
    @RV y ~ GaussianMeanVariance([2.0], mat(3.0))
    @RV z ~ GaussianMeanVariance([2.0], mat(3.0))
    n = Nonlinear{Sampling}(z, x, y, g=g, dims=[(1,), (1,), (1,)])

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    
    @test occursin("ruleSPNonlinearSInGX(g, 2, messages[3], messages[2], messages[1], dims=(1,))", code)
end

end # module