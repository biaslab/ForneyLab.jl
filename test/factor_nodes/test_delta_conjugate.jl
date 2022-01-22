module DeltaConjugateTest

using Test
using Random
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, logPdf, unsafeMean, unsafeVar, Conjugate, requiresBreaker, breakerParameters, renderCVI, naturalParams, ForgetDelayDescent
using ForneyLab: SPDeltaCInGX, SPDeltaCOutNM, SPDeltaCIn1MN, SPDeltaCOutNMX, SPDeltaCInMX, MDeltaCInMGX

f(x) = x
g(x) = x^2
h(x, y) = x + y

@testset "renderCVI" begin
    log_μ_bw(s) = -s'*s
    opt = ForgetDelayDescent(1.0, 1.0)
    
    # Multivariate Gaussian
    μ_fw = Message(Multivariate, GaussianMeanVariance, m=[1.0, 2.0], v=[2.0 1.0; 1.0 3.0])
    η = naturalParams(μ_fw.dist)
    λ = renderCVI(log_μ_bw, 1, opt, η, μ_fw)
    @test length(λ) == 6

    # Univariate Gaussian
    μ_fw = Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    η = naturalParams(μ_fw.dist)
    λ = renderCVI(log_μ_bw, 1, opt, η, μ_fw)
    @test length(λ) == 2

    # Gamma
    μ_fw = Message(Univariate, Gamma, a=1.0, b=2.0)
    η = naturalParams(μ_fw.dist)
    λ = renderCVI(log_μ_bw, 1, opt, η, μ_fw)
    @test length(λ) == 2

    # Categorical
    μ_fw = Message(Univariate, Categorical, p=[0.2, 0.3, 0.5])
    η = naturalParams(μ_fw.dist)
    λ = renderCVI(log_μ_bw, 1, opt, η, μ_fw)
    @test length(λ) == 3
end

@testset "requiresBreaker and breakerParameters" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Delta{Conjugate}(y, x, g=g)

    @test requiresBreaker(nd.i[:out]) # Single-input Delta{Conjugate} requires breaker
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    y = Variable()
    GaussianMeanVariance(z, 0.0, 1.0)
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Delta{Conjugate}(y, x, z, g=h)

    @test requiresBreaker(nd.i[:out]) # Multi-input Delta{Conjugate} requires breaker
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())
end


#-------------
# Update rules
#-------------

@testset "SPDeltaCOutNM" begin
    @test SPDeltaCOutNM <: SumProductRule{Delta{Conjugate}}
    @test outboundType(SPDeltaCOutNM) == Message{SampleList}
    @test isApplicable(SPDeltaCOutNM, [Nothing, Message{Gaussian}])
    @test isApplicable(SPDeltaCOutNM, [Nothing, Message{Bernoulli}])

    msg = ruleSPDeltaCOutNM(f, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny), n_samples=1)
    @test isapprox(msg.dist.params[:s][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
    msg = ruleSPDeltaCOutNM(f, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(tiny)), n_samples=1)
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "SPDeltaCIn1MN" begin
    @test SPDeltaCIn1MN <: SumProductRule{Delta{Conjugate}}
    @test outboundType(SPDeltaCIn1MN) == Message{FactorNode}
    @test isApplicable(SPDeltaCIn1MN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])

    msg = ruleSPDeltaCIn1MN(f, Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=0.1), n_samples=0)
    @test typeof(msg) == Message{GaussianWeightedMeanPrecision, Univariate}
end

@testset "SPDeltaCInGX" begin
    @test SPDeltaCInGX <: SumProductRule{Delta{Conjugate}}
    @test outboundType(SPDeltaCInGX) == Message{GaussianWeightedMeanPrecision}
    @test !isApplicable(SPDeltaCInGX, [Nothing, Message{Gamma}])
    @test isApplicable(SPDeltaCInGX, [Message{Gaussian}, Message{Gaussian}, Nothing])

    msg_out = Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    msg_in2 = Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0)

    res = ruleSPDeltaCInGX(h, 1, msg_out, msg_in1, msg_in2, n_samples=1000)

    @test isapprox(mean(res.dist), 0.9998084187686186, atol=0.1)
    @test isapprox(var(res.dist), 1.9999999999999982, atol=0.1)
end

@testset "SPDeltaCOutNMX" begin
    @test SPDeltaCOutNMX <: SumProductRule{Delta{Conjugate}}
    @test outboundType(SPDeltaCOutNMX) == Message{SampleList}
    @test !isApplicable(SPDeltaCOutNMX, [Nothing, Message{Gaussian}])
    @test isApplicable(SPDeltaCOutNMX, [Nothing, Message{Gaussian}, Message{Gamma}])
    @test isApplicable(SPDeltaCOutNMX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    msg = ruleSPDeltaCOutNMX(h, nothing, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0), n_samples=1)
    @test typeof(msg) == Message{SampleList, Univariate}
end

@testset "SPDeltaCInMX" begin
    @test SPDeltaCInMX <: SumProductRule{Delta{Conjugate}}
    @test outboundType(SPDeltaCInMX) == Message{FactorNode}
    @test !isApplicable(SPDeltaCInMX, [Message, Message{Gaussian}])
    @test isApplicable(SPDeltaCInMX, [Message, Nothing, Message{Gaussian}, Message{Gamma}])
    @test isApplicable(SPDeltaCInMX, [Message, Nothing, Message{PointMass}])
    @test isApplicable(SPDeltaCInMX, [Message, Message{PointMass}, Nothing])

    msg = ruleSPDeltaCInMX(h, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0))
    @test typeof(msg) == Message{GaussianWeightedMeanPrecision, Univariate}
end

@testset "MDeltaCInMGX" begin
    @test MDeltaCInMGX <: MarginalRule{Delta{Conjugate}}
    @test isApplicable(MDeltaCInMGX, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(MDeltaCInMGX, [Nothing, Message{Gaussian}])
    @test !isApplicable(MDeltaCInMGX, [Nothing, Message{Gaussian}, Message{Gamma}])

    dist = ruleMDeltaCInMGX(h, Message(Univariate, GaussianMeanVariance, m=3.0, v=0.1), Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0))
    @test typeof(dist) == Distribution{Multivariate, GaussianWeightedMeanPrecision}
end

end # module