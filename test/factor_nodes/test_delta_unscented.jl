module DeltaUnscentedTest

using Test
using LinearAlgebra
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights, Unscented, requiresBreaker, breakerParameters
using ForneyLab: SPDeltaUTOutNG, SPDeltaUTIn1GG, SPDeltaUTOutNGX, SPDeltaUTInGX, MDeltaUTInGX
using ForneyLab: unscentedStatistics, smoothRTS, smoothRTSMessage, collectStatistics, marginalizeGaussianMV, concatenateGaussianMV, split, isMultiIn

f(x) = x

g(x::Float64) = x^2 - 5.0
g(x::Vector{Float64}) = x.^2 .- 5.0
g_inv(y::Float64) = sqrt(y + 5.0)
g_inv(y::Vector{Float64}) = sqrt.(y .+ 5.0)

h(x::Float64, y::Float64) = x^2 - y
h(x::Vector{Float64}, y::Vector{Float64}) = x.^2 .- y
h(x::Float64, y::Vector{Float64}) = x^2 .- y
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

    # Mixed variate inbounds
    (m_tilde, V_tilde, C_tilde) = unscentedStatistics([m_x, [m_y, m_y]], [v_x, Diagonal([v_y, v_y])], h)
    @test m_tilde == [-1.0, -1.0]
    @test V_tilde == [5.000002000015229 2.000002000015229; 2.000002000015229 5.000002000015229]
    @test C_tilde == [0.0 0.0; -3.0000000000001137 0.0; 0.0 -3.0000000000001137]
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
    @test collectStatistics(Message(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0))) == ([[0.0], [2.0]], [mat(1.0), mat(3.0)])
    @test collectStatistics(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0))) == ([0.0, [2.0]], [1.0, mat(3.0)])
end

@testset "marginalizeGaussianMV" begin
    @test marginalizeGaussianMV([0.0, 1.0], [2.0 0.5; 0.5 3.0], [(), ()], 1) == (0.0, 2.0) # Univariate
    @test marginalizeGaussianMV([0.0, 1.0, 2.0], [2.0 0.0 0.5; 0.0 3.0 0.0; 0.5 0.0 4.0], [(), (2,)], 1) == (0.0, 2.0) # Univariate
    @test marginalizeGaussianMV([0.0, 1.0, 2.0], [2.0 0.0 0.5; 0.0 3.0 0.0; 0.5 0.0 4.0], [(), (2,)], 2) == ([1.0, 2.0], [3.0 0.0; 0.0 4.0]) # Multivariate
end

@testset "concatenateGaussianMV" begin
    @test concatenateGaussianMV([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) == ([1.0, 2.0, 3.0], Diagonal([4.0, 5.0, 6.0]), [(), (), ()])
    @test concatenateGaussianMV([[1.0], [2.0, 3.0]], [mat(4.0), Diagonal([5.0, 6.0])]) == ([1.0, 2.0, 3.0], [4.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 6.0], [(1,), (2,)])
    @test concatenateGaussianMV([1.0, [2.0, 3.0]], [4.0, Diagonal([5.0, 6.0])]) == ([1.0, 2.0, 3.0], [4.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 6.0], [(), (2,)])
end

@testset "split" begin
    @test split([1.0, 2.0], [(), ()]) == [1.0, 2.0]
    @test split([1.0, 2.0, 3.0, 4.0], [(2,), (2,)]) == [[1.0, 2.0], [3.0, 4.0]]
    @test split([1.0, 2.0, 3.0], [(), (2,)]) == [1.0, [2.0, 3.0]]
    @test split([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [(2,), (2,2)]) == [[1.0, 2.0], [3.0 5.0; 4.0 6.0]]
end

@testset "requiresBreaker and breakerParameters" begin
    # Without given inverse
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Delta{Unscented}(y, x, g=g)
    
    @test requiresBreaker(nd.i[:out])
    @test_throws Exception breakerParameters(nd.i[:out].partner)
    @test breakerParameters(nd.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())

    # With given inverse
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    nd = GaussianMeanVariance(x, 0.0, 1.0)
    Delta{Unscented}(y, x, g=g, g_inv=g)
    
    @test !requiresBreaker(nd.i[:out])
    @test_throws Exception breakerParameters(nd.i[:out])
end

@testset "isMultiIn" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    GaussianMeanVariance(x, 0.0, 1.0)
    nd = Delta{Unscented}(y, x, g=g)
    Clamp(y, 1.0)
    @test !isMultiIn(nd)

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    z = Variable()
    GaussianMeanVariance(x, 0.0, 1.0)
    GaussianMeanVariance(z, 0.0, 1.0)
    nd = Delta{Unscented}(y, x, z, g=g)
    Clamp(y, 1.0)
    @test isMultiIn(nd)

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    z = Variable()
    GaussianMeanVariance(x, 0.0, 1.0)
    Clamp(z, 1.0)
    nd = Delta{Unscented}(y, x, z, g=g)
    Clamp(y, 1.0)
    @test !isMultiIn(nd)
end


#-------------
# Update rules
#-------------

@testset "SPDeltaUTOutNG" begin
    @test SPDeltaUTOutNG <: SumProductRule{Delta{Unscented}}
    @test outboundType(SPDeltaUTOutNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPDeltaUTOutNG, [Nothing, Message{GaussianMeanVariance}]) 

    @test ruleSPDeltaUTOutNG(g, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0)) == Message(Univariate, GaussianMeanVariance, m=2.0000000001164153, v=66.00000000093132)
    @test ruleSPDeltaUTOutNG(g, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.0, v=66.0)
    @test ruleSPDeltaUTOutNG(g, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0000000001164153], v=mat(66.00000000093132))
    @test ruleSPDeltaUTOutNG(g, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(66.0))
end

@testset "SPDeltaUTOutNGX" begin
    @test SPDeltaUTOutNGX <: SumProductRule{Delta{Unscented}}
    @test outboundType(SPDeltaUTOutNGX) == Message{GaussianMeanVariance}
    @test !isApplicable(SPDeltaUTOutNGX, [Nothing, Message{Gaussian}]) 
    @test isApplicable(SPDeltaUTOutNGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test !isApplicable(SPDeltaUTOutNGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    @test ruleSPDeltaUTOutNGX(h, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=1.9999999997671694, v=67.00000899657607)
    @test ruleSPDeltaUTOutNGX(h, nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.9999999997671694], v=mat(67.00000899657607))
end

@testset "SPDeltaUTIn1GG" begin
    @test SPDeltaUTIn1GG <: SumProductRule{Delta{Unscented}}
    @test outboundType(SPDeltaUTIn1GG) == Message{GaussianMeanVariance}
    @test isApplicable(SPDeltaUTIn1GG, [Message{Gaussian}, Nothing]) 

    # Without given inverse
    @test ruleSPDeltaUTIn1GG(g, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.499999999868301, v=0.3125000002253504)
    @test ruleSPDeltaUTIn1GG(g, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.5, v=0.3125)
    @test ruleSPDeltaUTIn1GG(g, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.499999999868301], v=mat(0.31250000021807445))
    @test ruleSPDeltaUTIn1GG(g, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.5], v=mat(0.3125))

    # With given inverse
    @test ruleSPDeltaUTIn1GG(g, g_inv, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.6255032138433307, v=0.10796282966583703)
    @test ruleSPDeltaUTIn1GG(g, g_inv, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.6251028535207217, v=0.10968772603524787)
    @test ruleSPDeltaUTIn1GG(g, g_inv, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.6255032138433307], v=mat(0.10796282966583703))
    @test ruleSPDeltaUTIn1GG(g, g_inv, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.6251028535207217], v=mat(0.10968772603524787))
end

@testset "SPDeltaUTInGX" begin
    @test SPDeltaUTInGX <: SumProductRule{Delta{Unscented}}
    @test outboundType(SPDeltaUTInGX) == Message{Gaussian}
    @test !isApplicable(SPDeltaUTInGX, [Message{Gaussian}, Nothing]) 
    @test !isApplicable(SPDeltaUTInGX, [Nothing, Message{Gaussian}, Message{Gaussian}]) 
    @test isApplicable(SPDeltaUTInGX, [Message{Gaussian}, Nothing, Message{Gaussian}]) 

    # Without given inverse
    @test ruleSPDeltaUTInGX(h, 1, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=6.666665554127903, w=2.6666662216903676)
    @test ruleSPDeltaUTInGX(h, 1, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[6.666665554127903], w=mat(2.666666221690368))

    # With given inverse
    @test ruleSPDeltaUTInGX(h, h_inv_x, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.6187538476660848, v=0.14431487274498522)
    @test ruleSPDeltaUTInGX(h, h_inv_x, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.6187538476660848], v=mat(0.14431487274475785))
end

@testset "MDeltaUTInGX" begin
    @test MDeltaUTInGX <: MarginalRule{Delta{Unscented}}
    @test isApplicable(MDeltaUTInGX, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleMDeltaUTInGX(h, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=5.0, v=1.0)) == Distribution(Multivariate, GaussianMeanVariance, m=[2.3636363470614055, 4.9090909132334355], v=[0.2727273058237252 0.1818181735464949; 0.18181817354649488 0.9545454566127697])
    @test ruleMDeltaUTInGX(h, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(1.0))) == Distribution(Multivariate, GaussianMeanVariance, m=[2.3636363470609245, 4.909090913233555], v=[0.2727273058246874 0.18181817354625435; 0.18181817354625435 0.9545454566128299])
end


#------------
# Integration
#------------

@testset "Delta integration via UT with given inverse" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Delta{Unscented}(y, x, g=g, g_inv=g_inv)
    
    @test isa(n, Delta{Unscented})
    
    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTOutNG(g, nothing, messages[2])", code)
    @test !occursin("g_inv", code)

    # Backward; g_inv should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTIn1GG(g, g_inv, messages[2], nothing)", code)
end

@testset "Multi-argument integration via UT" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    @RV z ~ GaussianMeanVariance(5.0, 1.0)
    n = Delta{Unscented}(y, x, z, g=h, g_inv=[h_inv_x, nothing])
    
    # Forward; h_inv_x should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTOutNGX(h, nothing, messages[3], messages[1])", code)
    @test !occursin("h_inv_x", code)

    # Backward with given inverse; h_inv_x should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTInGX(h, h_inv_x, messages[3], nothing, messages[1])", code)

    # Backward without given inverse
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(z)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTInGX(h, 2, messages[3], messages[2], messages[1])", code)
    @test !occursin("h_inv_x", code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", code)
end

@testset "Delta integration via UT with given alpha" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Delta{Unscented}(y, x, g=g, alpha=1.0)
    
    # Forward; alpha should be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTOutNG(g, nothing, messages[1], alpha=1.0)", code)
end

@testset "Delta integration via UT without given inverse" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Delta{Unscented}(y, x, g=g)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTOutNG(g, nothing, messages[1])", code)
    @test !occursin("g_inv", code)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)
    @test occursin("ruleSPDeltaUTIn1GG(g, messages[2], messages[1])", code)
    @test !occursin("g_inv", code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", code)
end

end #module