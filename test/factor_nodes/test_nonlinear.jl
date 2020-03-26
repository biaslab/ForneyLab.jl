module NonlinearTest

using Test
using Random
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights, prod!, logPdf, unsafeMean, unsafeVar, ProbabilityDistribution, Unscented, ImportanceSampling
using ForneyLab: SPNonlinearUTOutNG, SPNonlinearUTIn1GG, SPNonlinearUTOutNGX, SPNonlinearUTInGX, SPNonlinearISIn1MN, SPNonlinearISOutNG, MNonlinearUTNGX

Random.seed!(1234)


#--------
# Helpers
#--------

@testset "sigmaPointsAndWeights" begin
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(0.0, 1.0, alpha=1e-3)
    @test sigma_points == [0.0, 0.0010000000000143778, -0.0010000000000143778]
    @test weights_m == [-999998.9999712444, 499999.9999856222, 499999.9999856222]
    @test weights_c == [-999995.9999722444, 499999.9999856222, 499999.9999856222]

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights([0.0], mat(1.0), alpha=1e-3)
    @test sigma_points == [[0.0], [0.0010000000000143778], [-0.0010000000000143778]]
    @test weights_m == [-999998.9999712444, 499999.9999856222, 499999.9999856222]
    @test weights_c == [-999995.9999722444, 499999.9999856222, 499999.9999856222]
end

@testset "unscentedStatistics" begin
    @test true == false
end

@testset "smoothRTS" begin
    @test true == false
end

@testset "collectStatistics" begin
    @test true == false
end

@testset "slice" begin
    @test true == false
end

@testset "pack" begin
    @test true == false
end

@testset "unpack" begin
    @test true == false
end


#-------------
# Update rules
#-------------

f(x) = x

g(x::Float64) = x^2 - 5.0
g(x::Vector{Float64}) = x.^2 .- 5.0
g_inv(y::Float64) = sqrt(y + 5.0)
g_inv(y::Vector{Float64}) = sqrt.(y .+ 5.0)

h(x::Float64, y::Float64) = x^2 - y
h(x::Vector{Float64}, y::Vector{Float64}) = x.^2 .- y
h_inv_x(z::Float64, y::Float64) = sqrt(z + y)
h_inv_x(z::Vector{Float64}, y::Vector{Float64}) = sqrt.(z .+ y)

@testset "SPNonlinearUTOutNG" begin
    @test SPNonlinearUTOutNG <: SumProductRule{Nonlinear{Unscented}}
    @test outboundType(SPNonlinearUTOutNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearUTOutNG, [Nothing, Message{Gaussian}]) 

    @test ruleSPNonlinearUTOutNG(g, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0)) == Message(Univariate, GaussianMeanVariance, m=2.0000000001164153, v=66.00000000093132)
    @test ruleSPNonlinearUTOutNG(g, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.0, v=66.0)
    @test ruleSPNonlinearUTOutNG(g, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0000000001164153], v=mat(66.00000000093132))
    @test ruleSPNonlinearUTOutNG(g, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(66.0))
end

@testset "SPNonlinearUTOutNGX" begin
    @test SPNonlinearUTOutNGX <: SumProductRule{Nonlinear{Unscented}}
    @test outboundType(SPNonlinearUTOutNGX) == Message{GaussianMeanVariance}
    @test !isApplicable(SPNonlinearUTOutNGX, [Nothing, Message{Gaussian}]) 
    @test isApplicable(SPNonlinearUTOutNGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test !isApplicable(SPNonlinearUTOutNGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    @test ruleSPNonlinearUTOutNGX(h, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=1.9999999997671694, v=67.00000899797305)
    @test ruleSPNonlinearUTOutNGX(h, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.9999999997671694], v=mat(67.00000899657607))
end

@testset "SPNonlinearISOutNG" begin
    samples = 2.0 .+ randn(100000)
    p_dist = ProbabilityDistribution(Univariate, SampleList, s=samples)

    @test SPNonlinearISOutNG <: SumProductRule{Nonlinear{ImportanceSampling}}
    @test outboundType(SPNonlinearISOutNG) == Message{SampleList}
    @test isApplicable(SPNonlinearISOutNG, [Nothing, Message{Gaussian}])

    @test abs(unsafeMean(ruleSPNonlinearISOutNG(f, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)).dist) - unsafeMean(p_dist)) < 0.2
    @test abs(unsafeVar(ruleSPNonlinearISOutNG(f, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)).dist) - unsafeVar(p_dist)) < 0.2
end

@testset "SPNonlinearUTIn1GG" begin
    @test SPNonlinearUTIn1GG <: SumProductRule{Nonlinear{Unscented}}
    @test outboundType(SPNonlinearUTIn1GG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearUTIn1GG, [Message{Gaussian}, Nothing]) 

    # Without given inverse
    @test ruleSPNonlinearUTIn1GG(g, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.499999999868301, v=0.3125000002253504)
    @test ruleSPNonlinearUTIn1GG(g, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.5, v=0.3125)
    @test ruleSPNonlinearUTIn1GG(g, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.499999999868301], v=mat(0.31250000021807445))
    @test ruleSPNonlinearUTIn1GG(g, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.5], v=mat(0.3125))

    # With given inverse
    @test ruleSPNonlinearUTIn1GG(g, g_inv, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.6255032138433307, v=0.10796282966583703)
    @test ruleSPNonlinearUTIn1GG(g, g_inv, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.6251028535207217, v=0.10968772603524787)
    @test ruleSPNonlinearUTIn1GG(g, g_inv, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.6255032138433307], v=mat(0.10796282966583703))
    @test ruleSPNonlinearUTIn1GG(g, g_inv, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.6251028535207217], v=mat(0.10968772603524787))
end

@testset "SPNonlinearUTInGX" begin
    @test SPNonlinearUTInGX <: SumProductRule{Nonlinear{Unscented}}
    @test outboundType(SPNonlinearUTInGX) == Message{GaussianMeanVariance}
    @test !isApplicable(SPNonlinearUTInGX, [Message{Gaussian}, Nothing]) 
    @test !isApplicable(SPNonlinearUTInGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test isApplicable(SPNonlinearUTInGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    # Without given inverse
    @test ruleSPNonlinearUTInGX(h, 1, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.4705882352587847, v=0.21799313500892237)
    @test ruleSPNonlinearUTInGX(h, 1, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.4705882352587847], v=mat(0.21799313501053375))

    # With given inverse
    @test ruleSPNonlinearUTInGX(h, h_inv_x, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.6187538476660848, v=0.14431487274498522)
    @test ruleSPNonlinearUTInGX(h, h_inv_x, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.6187538476660848], v=mat(0.14431487274475785))
end

@testset "SPNonlinearISIn1MN" begin
    @test SPNonlinearISIn1MN <: SumProductRule{Nonlinear{ImportanceSampling}}
    @test outboundType(SPNonlinearISIn1MN) == Message{Function}
    @test isApplicable(SPNonlinearISIn1MN, [Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing])
    
    log_pdf(x) = ruleSPNonlinearISIn1MN(f, Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), nothing).dist.params[:log_pdf](x)
    @test log_pdf(1.5) == logPdf(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=1.0), 1.5)
end

@testset "MNonlinearUTNGX" begin
    @test MNonlinearUTNGX <: MarginalRule{Nonlinear{Unscented}}
    @test isApplicable(MNonlinearUTNGX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleMNonlinearUTNGX(h, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.3636363470614056, -0.09090908676656449], v=[0.2727273058237252 0.1818181735464949; 0.18181817354649488 0.9545454566127697])
    @test ruleMNonlinearUTNGX(h, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.36363634706092446, -0.09090908676644421], v=[0.2727273058246874 0.18181817354625435; 0.18181817354625435 0.9545454566128299])
end


#------------
# Integration
#------------

@testset "Nonlinear integration via UT with given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g=g, g_inv=g_inv)
    
    @test isa(n, Nonlinear{Unscented})
    
    # Forward; g_inv should not be present in call
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTOutNG(g, nothing, messages[2])", algo_code)
    @test !occursin("g_inv", algo_code)

    # Backward; g_inv should be present in call
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTIn1GG(g, g_inv, messages[2], nothing)", algo_code)
end

@testset "Multi-argument nonlinear integration via UT" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    @RV z ~ GaussianMeanVariance(5.0, 1.0)
    n = Nonlinear(y, x, z, g=h, g_inv=[h_inv_x, nothing])
    
    # Forward; h_inv_x should not be present in call
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTOutNGX(h, nothing, messages[2], messages[3])", algo_code)
    @test !occursin("h_inv_x", algo_code)

    # Backward with given inverse; h_inv_x should be present in call
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTInGX(h, h_inv_x, messages[2], nothing, messages[3])", algo_code)

    # Backward without given inverse
    algo = sumProductAlgorithm(z)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTInGX(h, 2, messages[3], messages[2], messages[1])", algo_code)
    @test !occursin("h_inv_x", algo_code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", algo_code)
end

@testset "Nonlinear integration via UT with given alpha" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{Unscented}(y, x, g=g, alpha=1.0)
    
    # Forward; alpha should be present in call
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTOutNG(g, nothing, messages[2], alpha=1.0)", algo_code)
end

@testset "Nonlinear integration via UT without given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{Unscented}(y, x, g=g)

    # Forward; g_inv should not be present in call
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTOutNG(g, nothing, messages[2])", algo_code)
    @test !occursin("$(string(g_inv))", algo_code)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearUTIn1GG(g, messages[2], messages[1])", algo_code)
    @test !occursin("g_inv", algo_code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", algo_code)
end

@testset "Nonlinear integration via importance sampling" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear{ImportanceSampling}(y, x, g=g)

    # Forward; g_inv should not be present in call
    algo = sumProductAlgorithm(y)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearISOutNG(g, nothing, messages[2])", algo_code)
end

end #module