module NonlinearTest

using Test
using Random
using LinearAlgebra
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights, prod!, logPdf, unsafeMean, unsafeVar, Unscented
using ForneyLab: SPNonlinearUTOutNG, SPNonlinearUTIn1GG, SPNonlinearUTOutNGX, SPNonlinearUTInGX, MNonlinearUTInGX
using ForneyLab: unscentedStatistics, smoothRTS, smoothRTSMessage, collectStatistics, marginalizeGaussianMV, concatenateGaussianMV, split

Random.seed!(1234)

f(x) = x

g(x::Float64) = x^2 - 5.0
g(x::Vector{Float64}) = x.^2 .- 5.0
g_inv(y::Float64) = sqrt(y + 5.0)
g_inv(y::Vector{Float64}) = sqrt.(y .+ 5.0)

h(x::Float64, y::Float64) = x^2 - y
h(x::Vector{Float64}, y::Vector{Float64}) = x.^2 .- y
h_inv_x(z::Float64, y::Float64) = sqrt(z + y)
h_inv_x(z::Vector{Float64}, y::Vector{Float64}) = sqrt.(z .+ y)


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
    m_x = 0.0
    v_x = 1.0
    m_y = 2.0
    v_y = 3.0

    # Single univariate inbound
    (m_tilde, V_tilde, C_tilde) = unscentedStatistics(m_x, v_x, g)
    @test m_tilde == -4.0
    @test V_tilde == 1.9999999997671694
    @test C_tilde == 0.0

    # Multiple univariate inbounds
    (m_tilde, V_tilde, C_tilde) = unscentedStatistics([m_x, m_y], [v_x, v_y], h)
    @test m_tilde == -1.0000000000582077
    @test V_tilde == 5.0000009994837455
    @test C_tilde == [0.0, -2.9999999999442934]

    # Single multivariate inbound
    (m_tilde, V_tilde, C_tilde) = unscentedStatistics([m_x, m_y], Diagonal([v_x, v_y]), g)
    @test m_tilde == [-4.000000000465661, 1.9999999997380655]
    @test V_tilde == [2.000000997039024 5.999996995087713; 5.999996995087713 66.00000899611041]
    @test C_tilde == [0.0 0.0; 5.547917680814862e-11 12.000000000165528]

    # Multiple multivariate inbounds
    (m_tilde, V_tilde, C_tilde) = unscentedStatistics([[m_x, m_x], [m_y, m_y]], [Diagonal([v_x, v_x]), Diagonal([v_y, v_y])], h)
    @test m_tilde == [-1.000000000174623, -1.000000000174623]
    @test V_tilde == [5.000002998916898 1.9999989990901668; 1.999998999031959 5.000002998946002]
    @test C_tilde == [0.0 0.0; 0.0 0.0; -2.9999999999998863 0.0; 0.0 -2.9999999999998863]
end

@testset "smoothRTSMessage" begin
    @test smoothRTSMessage(-4.0, 2.0, 1.0, 3.0, 5.0, 4.0, 6.0) == (43.0, 195.0)
    @test smoothRTSMessage([-4.0], mat(2.0), mat(1.0), [3.0], mat(5.0), [4.0], mat(6.0)) == ([43.0], mat(195.0))
end

@testset "smoothRTS" begin
    @test smoothRTS(-4.0, 2.0, 1.0, 3.0, 5.0, 4.0, 6.0) == (4.0, 4.875)
    @test smoothRTS([-4.0], mat(2.0), mat(1.0), [3.0], mat(5.0), [4.0], mat(6.0)) == ([4.0], mat(4.875))
end

@testset "collectStatistics" begin
    @test collectStatistics(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0)) == ([0.0, 2.0], [1.0, 3.0])
    @test collectStatistics(Message(Univariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)), nothing, Message(Univariate, GaussianMeanVariance, m=[2.0], v=mat(3.0))) == ([[0.0], [2.0]], [mat(1.0), mat(3.0)])
end

@testset "marginalizeGaussianMV" begin
    @test marginalizeGaussianMV(Univariate, [0.0, 1.0], [2.0 0.5; 0.5 3.0], ones(Int64, 2), 1) == (0.0, 2.0)
    @test marginalizeGaussianMV(Multivariate, [0.0, 1.0, 2.0], [2.0 0.0 0.5; 0.0 3.0 0.0; 0.5 0.0 4.0], [1, 2], 2) == ([1.0, 2.0], [3.0 0.0; 0.0 4.0])
end

@testset "concatenateGaussianMV" begin
    @test concatenateGaussianMV([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) == ([1.0, 2.0, 3.0], Diagonal([4.0, 5.0, 6.0]), ones(Int64, 3))
    @test concatenateGaussianMV([[1.0], [2.0, 3.0]], [mat(4.0), Diagonal([5.0, 6.0])]) == ([1.0, 2.0, 3.0], [4.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 6.0], [1, 2])
end

@testset "split" begin
    @test split([1.0, 2.0, 3.0], [1, 2]) == [[1.0], [2.0, 3.0]]
end


#-------------
# Update rules
#-------------

@testset "SPNonlinearUTOutNG" begin
    @test SPNonlinearUTOutNG <: SumProductRule{Nonlinear{Unscented}}
    @test outboundType(SPNonlinearUTOutNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearUTOutNG, [Nothing, Message{GaussianMeanVariance}]) 

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
    @test outboundType(SPNonlinearUTInGX) == Message{Gaussian}
    @test !isApplicable(SPNonlinearUTInGX, [Message{Gaussian}, Nothing]) 
    @test !isApplicable(SPNonlinearUTInGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test isApplicable(SPNonlinearUTInGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    # Without given inverse
    @test ruleSPNonlinearUTInGX(h, 1, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=6.666665554160243, w=2.6666662217033044)
    @test ruleSPNonlinearUTInGX(h, 1, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[6.666665554127903], w=mat(2.666666221690368))

    # With given inverse
    @test ruleSPNonlinearUTInGX(h, h_inv_x, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.6187538476660848, v=0.14431487274498522)
    @test ruleSPNonlinearUTInGX(h, h_inv_x, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.6187538476660848], v=mat(0.14431487274475785))
end

@testset "MNonlinearUTInGX" begin
    @test MNonlinearUTInGX <: MarginalRule{Nonlinear{Unscented}}
    @test isApplicable(MNonlinearUTInGX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleMNonlinearUTInGX(h, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.3636363470614055, 4.9090909132334355], v=[0.2727273058237252 0.1818181735464949; 0.18181817354649488 0.9545454566127697])
    @test ruleMNonlinearUTInGX(h, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.3636363470609245, 4.909090913233555], v=[0.2727273058246874 0.18181817354625435; 0.18181817354625435 0.9545454566128299])
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

end #module