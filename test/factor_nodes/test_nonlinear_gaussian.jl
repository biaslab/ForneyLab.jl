module NonlinearGaussianTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, sigmaPointsAndWeights
using ForneyLab: SPNonlinearGaussianOutNGP, SPNonlinearGaussianIn1GGP, MNonlinearGaussianGGD

@testset "sigmaPointsAndWeights" begin
    dist = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(dist, alpha=1e-3)
    @test sigma_points == [0.0, 0.0010000000000143778, -0.0010000000000143778]
    @test weights_m == [-999998.9999712444, 499999.9999856222, 499999.9999856222]
    @test weights_c == [-999995.9999722444, 499999.9999856222, 499999.9999856222]

    dist = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(dist, alpha=1e-3)
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

@testset "SPNonlinearGaussianOutNGP" begin
    @test SPNonlinearGaussianOutNGP <: SumProductRule{NonlinearGaussian}
    @test outboundType(SPNonlinearGaussianOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearGaussianOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, PointMass, m=1.0), g) == Message(Univariate, GaussianMeanVariance, m=2.0000000001164153, v=67.00000000093132)
    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, PointMass, m=1.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.0, v=67.0)
    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.0000000001164153], v=mat(67.00000000093132))
    @test ruleSPNonlinearGaussianOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(67.0))
end

@testset "SPNonlinearGaussianIn1GGP" begin
    @test SPNonlinearGaussianIn1GGP <: SumProductRule{NonlinearGaussian}
    @test outboundType(SPNonlinearGaussianIn1GGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPNonlinearGaussianIn1GGP, [Message{Gaussian}, Nothing, Message{PointMass}]) 

    # Without given inverse
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, PointMass, m=1.0), g) == Message(Univariate, GaussianMeanVariance, m=2.499999999868301, v=0.3750000002288054)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), Message(Univariate, GaussianMeanVariance, m=2.0, v=1.0), Message(Univariate, PointMass, m=1.0), g, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.5, v=0.375)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g) == Message(Multivariate, GaussianMeanVariance, m=[2.499999999868301], v=mat(0.37500000022152946))
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(1.0)), g, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.5], v=mat(0.375))

    # With given inverse
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, PointMass, m=1.0), g, g_inv) == Message(Univariate, GaussianMeanVariance, m=2.6187538478989154, v=0.14431487180792146)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0), nothing, Message(Univariate, PointMass, m=1.0), g, g_inv, alpha=1.0) == Message(Univariate, GaussianMeanVariance, m=2.618033988749895, v=0.14743453366290887)
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(1.0)), g, g_inv) == Message(Multivariate, GaussianMeanVariance, m=[2.6187538478989154], v=mat(0.14431487180792146))
    @test ruleSPNonlinearGaussianIn1GGP(Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(1.0)), g, g_inv, alpha=1.0) == Message(Multivariate, GaussianMeanVariance, m=[2.618033988749895], v=mat(0.14743453366290887))
end

@testset "MNonlinearGaussianGGD" begin
    @test MNonlinearGaussianGGD <: MarginalRule{NonlinearGaussian}
    @test isApplicable(MNonlinearGaussianGGD, [Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution])

    @test ruleMNonlinearGaussianGGD(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), ProbabilityDistribution(Univariate, PointMass, m=1.0), g) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[-0.4444444446018769, 1.388888889228411], w=[1.1111111111424605 -0.22222222229108557; -0.22222222229108557 0.9444444445944998])
end

@testset "averageEnergy" begin
    @test averageEnergy(NonlinearGaussian, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)), ProbabilityDistribution(Univariate, PointMass, m=1.0), g) == 10.418938533204672
end


#------------
# Integration
#------------

@testset "NonlinearGaussian integration with given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = NonlinearGaussian(y, x, 1.0, g, g_inv=g_inv)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization()
    algo = sumProductAlgorithm(y, pfz)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearGaussianOutNGP(nothing, messages[2], Message(Univariate, PointMass, m=1.0), g)", algo_code)
    @test !occursin("g_inv", algo_code)

    # Backward; g_inv should be present in call
    pfz = PosteriorFactorization()
    algo = sumProductAlgorithm(x, pfz)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearGaussianIn1GGP(messages[2], nothing, Message(Univariate, PointMass, m=1.0), g, g_inv)", algo_code)
end

@testset "NonlinearGaussian integration with given alpha" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = NonlinearGaussian(y, x, 1.0, g, alpha=1.0)

    # Forward; alpha should be present in call
    pfz = PosteriorFactorization()
    algo = sumProductAlgorithm(y, pfz)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearGaussianOutNGP(nothing, messages[2], Message(Univariate, PointMass, m=1.0), g, alpha=1.0)", algo_code)
end

@testset "NonlinearGaussian integration without given inverse" begin
    FactorGraph()

    @RV x ~ GaussianMeanVariance(2.0, 1.0)
    @RV y ~ GaussianMeanVariance(2.0, 3.0)
    n = NonlinearGaussian(y, x, 1.0, g)

    # Forward; g_inv should not be present in call
    pfz = PosteriorFactorization()
    algo = sumProductAlgorithm(y, pfz)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearGaussianOutNGP(nothing, messages[2], Message(Univariate, PointMass, m=1.0), g)", algo_code)
    @test !occursin("$(string(g_inv))", algo_code)

    # Backward; g_inv should not be present in call, 
    # both messages should be required, and initialization should take place
    pfz =  PosteriorFactorization()
    algo = sumProductAlgorithm(x, pfz)
    algo_code = algorithmSourceCode(algo)
    @test occursin("ruleSPNonlinearGaussianIn1GGP(messages[2], messages[1], Message(Univariate, PointMass, m=1.0), g)", algo_code)
    @test !occursin("g_inv", algo_code)
    @test occursin("messages[1] = Message(vague(GaussianMeanVariance))", algo_code)
end

end # module