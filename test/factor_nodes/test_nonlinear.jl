module NonlinearTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights
import ForneyLab: SPNonlinearOutNG, SPNonlinearIn1GG

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

@testset "SPNonlinearOutNG" begin
    @test SPNonlinearOutNG <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearOutNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearOutNG, [Nothing, Message{Gaussian}]) 

    @test ruleSPNonlinearOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), g) == Message(Univariate, GaussianMeanVariance, m=2.0000000001164153, v=66.00000000093132)
    @test ruleSPNonlinearOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.0, v=66.0)
    @test ruleSPNonlinearOutNG(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.0000000001164153], v=mat(66.00000000093132))
    @test ruleSPNonlinearOutNG(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(66.0))
end

@testset "SPNonlinearIn1GG" begin
    @test SPNonlinearIn1GG <: SumProductRule{Nonlinear}
    @test outboundType(SPNonlinearIn1GG) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearIn1GG, [Message{Gaussian}, Nothing]) 

    # Without given inverse
    @test ruleSPNonlinearIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), g) == Message(Univariate, GaussianMeanVariance, m=2.499999999868301, v=0.3125000002253504)
    @test ruleSPNonlinearIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.5, v=0.3125)
    @test ruleSPNonlinearIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.499999999868301], v=mat(0.31250000021807445))
    @test ruleSPNonlinearIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.5], v=mat(0.3125))

    # With given inverse
    @test ruleSPNonlinearIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, g, g_inv) == Message(Univariate, GaussianMeanVariance, m=2.6255032138433307, v=0.10796282966583703)
    @test ruleSPNonlinearIn1GG(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, g, g_inv, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.6251028535207217, v=0.10968772603524787)
    @test ruleSPNonlinearIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, g, g_inv) == Message(Multivariate, GaussianMeanVariance, m=[2.6255032138433307], v=mat(0.10796282966583703))
    @test ruleSPNonlinearIn1GG(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, g, g_inv, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.6251028535207217], v=mat(0.10968772603524787))
end


#------------
# Integration
#------------

@testset "Nonlinear integration with given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g, g_inv=g_inv)

    # Forward; g_inv should not be present in call
    algo = sumProductAlgorithm(y)
    @test occursin("ruleSPNonlinearOutNG(nothing, messages[2], $(string(g))", algo)
    @test !occursin("g_inv", algo)

    # Backward; g_inv should be present in call
    algo = sumProductAlgorithm(x)
    @test occursin("ruleSPNonlinearIn1GG(messages[2], nothing, $(string(g)), $(string(g_inv)))", algo)
end

@testset "Nonlinear integration with given alpha" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g, alpha=1.0)

    # Forward; alpha should be present in call
    algo = sumProductAlgorithm(y)
    @test occursin("ruleSPNonlinearOutNG(nothing, messages[2], $(string(g)), alpha=1.0)", algo)
end

@testset "Nonlinear integration without given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = Nonlinear(y, x, g)

    # Forward; g_inv should not be present in call
    algo = sumProductAlgorithm(y)
    @test occursin("ruleSPNonlinearOutNG(nothing, messages[2], $(string(g)))", algo)
    @test !occursin("$(string(g_inv))", algo)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    algo = sumProductAlgorithm(x)
    @test occursin("ruleSPNonlinearIn1GG(messages[2], messages[1], $(string(g)))", algo)
    @test !occursin("g_inv", algo)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", algo)
end

end # module