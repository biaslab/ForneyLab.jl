module NonlinearTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights, prod!, logPdf, unsafeMean, unsafeVar, ProbabilityDistribution, SPNonlinearUTOutNG, SPNonlinearUTIn1GG, SPNonlinearPTInMN, SPNonlinearPTOutNG

@testset "sigmaPointsAndWeights" begin
    dist = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(dist, 1e-3)
    @test sigma_points == [0.0, 0.0010000000000143778, -0.0010000000000143778]
    @test weights_m == [-999998.9999712444, 499999.9999856222, 499999.9999856222]
    @test weights_c == [-999995.9999722444, 499999.9999856222, 499999.9999856222]

    dist = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(dist, 1e-3)
    @test sigma_points == [[0.0], [0.0010000000000143778], [-0.0010000000000143778]]
    @test weights_m == [-999998.9999712444, 499999.9999856222, 499999.9999856222]
    @test weights_c == [-999995.9999722444, 499999.9999856222, 499999.9999856222]
end


#-------------
# Update rules
#-------------

g(x::Float64) = x^2 - 5.0
g(x::Vector{Float64}) = x.^2 .- 5.0
g_inv(x::Float64) = sqrt(x + 5.0)
g_inv(x::Vector{Float64}) = sqrt.(x .+ 5.0)

@testset "SPNonlinearUTOutNG" begin
    @test SPNonlinearUTOutNG <: SumProductRule{NonlinearUT}
    @test outboundType(SPNonlinearUTOutNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearUTOutNG, [Nothing, Message{Gaussian}]) 

    @test ruleSPNonlinearUTOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), g) == Message(Univariate, GaussianMeanVariance, m=2.0000000001164153, v=66.00000000093132)
    @test ruleSPNonlinearUTOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.0, v=66.0)
    @test ruleSPNonlinearUTOutNG(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.0000000001164153], v=mat(66.00000000093132))
    @test ruleSPNonlinearUTOutNG(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(66.0))
end

@testset "SPNonlinearUTIn1GG" begin
    @test SPNonlinearUTIn1GG <: SumProductRule{NonlinearUT}
    @test outboundType(SPNonlinearUTIn1GG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearUTIn1GG, [Message{Gaussian}, Nothing]) 

    # Without given inverse
    @test ruleSPNonlinearUTIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), g) == Message(Univariate, GaussianMeanVariance, m=2.499999999868301, v=0.3125000002253504)
    @test ruleSPNonlinearUTIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.5, v=0.3125)
    @test ruleSPNonlinearUTIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.499999999868301], v=mat(0.31250000021807445))
    @test ruleSPNonlinearUTIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.5], v=mat(0.3125))

    # With given inverse
    @test ruleSPNonlinearUTIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, g, g_inv) == Message(Univariate, GaussianMeanVariance, m=2.6255032138433307, v=0.10796282966583703)
    @test ruleSPNonlinearUTIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, g, g_inv, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.6251028535207217, v=0.10968772603524787)
    @test ruleSPNonlinearUTIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, g, g_inv) == Message(Multivariate, GaussianMeanVariance, m=[2.6255032138433307], v=mat(0.10796282966583703))
    @test ruleSPNonlinearUTIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, g, g_inv, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.6251028535207217], v=mat(0.10968772603524787))
end


#------------
# Integration
#------------

@testset "Nonlinear integration via UT with given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g)
    applyUnscentedTransform(y, g_inv=g_inv)
    
    # Forward; g_inv should not be present in call
    algo = InferenceAlgorithm()
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearOutNG(nothing, messages[2], g)", algo_code)
    @test !occursin("g_inv", algo_code)

    # Backward; g_inv should be present in call
    algo = InferenceAlgorithm()
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearIn1GG(messages[2], nothing, g, g_inv)", algo_code)
end

@testset "Nonlinear integration via UT with given alpha" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g)
    applyUnscentedTransform(y, alpha=1.0)
    
    # Forward; alpha should be present in call
    algo = InferenceAlgorithm()
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearOutNG(nothing, messages[2], g, alpha=1.0)", algo_code)
end

@testset "Nonlinear integration via UT without given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g)
    applyUnscentedTransform(y)

    # Forward; g_inv should not be present in call
    algo = InferenceAlgorithm()
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearOutNG(nothing, messages[2], g)", algo_code)
    @test !occursin("$(string(g_inv))", algo_code)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    algo =  InferenceAlgorithm()
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearIn1GG(messages[2], messages[1], g)", algo_code)
    @test !occursin("g_inv", algo_code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", algo_code)
end

@testset "prod!" begin
    f_dummy(x) = x
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0))).params[:m]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0)).params[:m]) < 0.1
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0))).params[:v]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0)).params[:v]) < 0.1
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.2, v=1.0)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.6, v=0.5))).params[:m]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=1.2, v=1.0),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.6, v=0.5)).params[:m]) < 0.1
    @test abs((convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=6.5, v=4.1)
            *ProbabilityDistribution(Univariate, GaussianMeanVariance, m=12.0, v=3.0))).params[:v]
            -  (ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=6.5, v=4.1),nothing,f_dummy).dist*ProbabilityDistribution(Univariate, GaussianMeanVariance, m=12.0, v=3.0)).params[:v]) < 0.1
end

#-------------
# Update rules
#-------------

@testset "SPNonlinearPTInMN" begin
    f_dummy(x) = x
    @test SPNonlinearPTInMN <: SumProductRule{NonlinearPT}
    @test outboundType(SPNonlinearPTInMN) == Message{Function}
    @test isApplicable(SPNonlinearPTInMN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])
    f(x) = ruleSPNonlinearPTInMN(Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),nothing,f_dummy).dist.params[:log_pdf](x)
    @test f(1.5) == logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), 1.5)
end

 @testset "SPNonlinearPTOutNG" begin
     f_dummy(x) = x
     samples = 2.0 .+ randn(100000)
     p_dist = ProbabilityDistribution(Univariate, SampleList, s=samples)
     @test SPNonlinearPTOutNG <: SumProductRule{NonlinearPT}
     @test outboundType(SPNonlinearPTOutNG) == Message{SampleList}
     @test isApplicable(SPNonlinearPTOutNG, [Nothing, Message{Gaussian}])
     @test abs(unsafeMean(ruleSPNonlinearPTOutNG(nothing,Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),f_dummy).dist) - unsafeMean(p_dist)) < 0.2
     @test abs(unsafeVar(ruleSPNonlinearPTOutNG(nothing,Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0),f_dummy).dist) - unsafeVar(p_dist)) < 0.2
 end

end #module