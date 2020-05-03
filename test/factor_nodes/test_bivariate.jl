module BivariateTest

using Test
using Random
# using LinearAlgebra
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights, prod!, logPdf, unsafeMean, unsafeVar, ProbabilityDistribution, Unscented, Sampling
using ForneyLab: SPBivariateSIn1MNG, ruleSPBivariateSIn1MNG, SPBivariateSIn2MGN, ruleSPBivariateSIn2MGN, SPBivariateSOutNGG, ruleSPBivariateSOutNGG, mergeInputs, gradientOptimization, decomposePosteriorParameters
using ForwardDiff

Random.seed!(1234)

g(x) = x^2

@testset "mergeInputs" begin
    dist1 = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    dist2 = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    dist3 = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0, 3.0], v=[1.0 0.0; 0.0 1.0])
    dist4 = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[4.0, 5.0], v=[2.0 0.0; 0.0 3.0])

    # Univariate + Univariate
    m_concat, v_concat, dim1, dim2 = mergeInputs(dist1, dist2)
    @test m_concat == [0.0, 0.0]
    @test v_concat == [1.0 0.0; 0.0 1.0]
    @test dim1 == 1
    @test dim2 == 1
    
    # Univariate + Multivariate
    m_concat, v_concat, dim1, dim2 = mergeInputs(dist1, dist3)
    @test m_concat == [0.0, 2.0, 3.0]
    @test v_concat == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test dim1 == 1
    @test dim2 == 2
    
    # Multivariate + Multivariate
    m_concat, v_concat, dim1, dim2 = mergeInputs(dist3, dist4)
    @test m_concat == [2.0, 3.0, 4.0, 5.0]
    @test v_concat == [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 2.0 0.0; 0.0 0.0 0.0 3.0]
    @test dim1 == 2
    @test dim2 == 2
    
end

@testset "decomposePosteriorDistributions" begin
    dist1 = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    dist2 = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    dist3 = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0, 3.0], v=[1.0 0.0; 0.0 1.0])
    dist4 = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[4.0, 5.0], v=[2.0 0.0; 0.0 3.0])

    # Univariate + Univariate
    m_concat, v_concat, dim1, dim2 = mergeInputs(dist1, dist2)
    m1, v1, m2, v2 = decomposePosteriorParameters(dist1, dist2, m_concat, v_concat)
    @test (m1, v1, m2, v2) == (0.0, 1.0, 0.0, 1.0)
    
    # Univariate + Multivariate
    m_concat, v_concat, dim1, dim2 = mergeInputs(dist1, dist3)
    m1, v1, m2, v2 = decomposePosteriorParameters(dist1, dist3, m_concat, v_concat)
    @test (m1, v1, m2, v2) == (0.0, 1.0, [2.0, 3.0], [1.0 0.0; 0.0 1.0])
    
    # Multivariate + Multivariate
    m_concat, v_concat, dim1, dim2 = mergeInputs(dist3, dist4)
    m1, v1, m2, v2 = decomposePosteriorParameters(dist3, dist4, m_concat, v_concat)
    @test (m1, v1, m2, v2) == ([2.0, 3.0], [1.0 0.0; 0.0 1.0], [4.0, 5.0], [2.0 0.0; 0.0 3.0])
    
end


@testset "gradientOptimization" begin

f(x) = -x[1]^2-x[2]^2
grad(x) = ForwardDiff.gradient(f, x)
res = gradientOptimization(f, grad, [4.0, 1.0], 0.001)

# Optimum is [0.0, 0.0]
@test isapprox(res, [0.0, 0.0], atol=1e-16)

end

@testset "ruleSPBivariateSOutNGG" begin
    @test SPBivariateSOutNGG <: SumProductRule{Bivariate{Sampling}}
    @test outboundType(SPBivariateSOutNGG) == Message{SampleList}
    @test !isApplicable(SPBivariateSOutNGG, [Nothing, Message{Poisson}]) 
    @test isApplicable(SPBivariateSOutNGG, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    
    h(x,y) = x+y

    msg_in1 = Message(GaussianMeanVariance, m=0.0, v=1.0)
    msg_in2 = Message(GaussianMeanVariance, m=2.0, v=1.0)
    
    res = ruleSPBivariateSOutNGG(nothing, msg_in1, msg_in2, h, Dict(), 100000)
    @test isapprox(mean(res.dist), 2.0, atol=0.1)
    @test isapprox(var(res.dist), 2.0, atol=0.1)
end



@testset "ruleSPBivariateSIn1MNG" begin
    @test SPBivariateSIn1MNG <: SumProductRule{Bivariate{Sampling}}
    @test outboundType(SPBivariateSIn1MNG) == Message{Gaussian}
    @test !isApplicable(SPBivariateSIn1MNG, [Nothing, Message{Gaussian}]) 
    @test isApplicable(SPBivariateSIn1MNG, [Message{FactorFunction}, Nothing, Message{Gaussian}]) 
    
    msg_out = Message(vague(GaussianMeanVariance))
    msg_in1 = Message(vague(GaussianMeanVariance))
    msg_in2 = Message(vague(GaussianMeanVariance))
    status = Dict{Symbol, Any}()
    status[:updated] = true
    status[:message] = Message(GaussianMeanVariance,m=2.0, v=4.0)
    
    res = ruleSPBivariateSIn1MNG(msg_out, msg_in1, msg_in2, g, status, 1000)

    @test mean(res.dist) == 2.0
    @test var(res.dist) == 4.0
    @test !(status[:updated])
end

@testset "ruleSPBivariateSIn2MGN" begin
    @test SPBivariateSIn2MGN <: SumProductRule{Bivariate{Sampling}}
    @test outboundType(SPBivariateSIn2MGN) == Message{Gaussian}
    @test !isApplicable(SPBivariateSIn2MGN, [Nothing, Message{Gamma}]) 
    @test isApplicable(SPBivariateSIn2MGN, [Message{FactorFunction}, Message{Gaussian}, Nothing]) 
    
    msg_out = Message(vague(GaussianMeanVariance))
    msg_in1 = Message(vague(GaussianMeanVariance))
    msg_in2 = Message(vague(GaussianMeanVariance))
    status = Dict{Symbol, Any}()
    status[:updated] = true
    status[:message] = Message(GaussianMeanVariance,m=2.0, v=4.0)
    
    res = ruleSPBivariateSIn2MGN(msg_out, msg_in1, msg_in2, g, status, 1000)

    @test mean(res.dist) == 2.0
    @test var(res.dist) == 4.0
    @test !(status[:updated])
end

@testset "Bivariate integration via sampling" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    @RV z ~ GaussianMeanVariance(2.0, 3.0)
    n = Bivariate{Sampling}(z, x, y, g)

    # Forward; g_inv should not be present in call
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    
    @test occursin("ruleSPBivariateSIn2MGN(messages[3], messages[2], messages[1], g, currentGraph().nodes[:bivariate_1].status, 1000)", algo_code)
end

end