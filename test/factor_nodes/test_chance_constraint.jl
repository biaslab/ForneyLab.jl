module ChanceConstraintTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPChanceConstraintOutG, standardGaussianPdf, standardGaussianCdf, truncatedGaussianMoments


#--------
# Helpers
#--------

@testset "standardGaussianPdf" begin
    @test standardGaussianPdf(-Inf) == 0.0
    @test standardGaussianPdf(0.0) == 0.3989422804014327
    @test standardGaussianPdf(Inf) == 0.0
end

@testset "standardGaussianPdf" begin
    @test standardGaussianCdf(-Inf) == 0.0
    @test standardGaussianCdf(0.0) == 0.5
    @test standardGaussianCdf(Inf) == 1.0
end

@testset "truncatedGaussianMoments" begin
    @test truncatedGaussianMoments(0.0, 1.0, 0.0, Inf) == (0.5, 0.7978845608028654, 0.3633802276324186)
    @test truncatedGaussianMoments(0.0, 1.0, -Inf, 0.0) == (0.5, -0.7978845608028654, 0.3633802276324186)
    @test truncatedGaussianMoments(0.0, 1.0, -Inf, Inf) == (1.0, 0.0, 1.0)
    @test truncatedGaussianMoments(0.0, 1.0, 0.0, 0.0) == (0.0, 0.0, 0.0) # Improper truncation
end


#-------------
# Update rules
#-------------

@testset "SPChanceConstraintOutG" begin
    @test SPChanceConstraintOutG <: SumProductRule{ChanceConstraint}
    @test outboundType(SPChanceConstraintOutG) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPChanceConstraintOutG, [Nothing]) 

    @test ruleSPChanceConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), (0.0, Inf), 0.05) == Message(Univariate, GaussianWeightedMeanPrecision, xi=1.4826342923288633, w=1.064673910049474)
    @test ruleSPChanceConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), (-Inf, Inf), 0.05) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
end


#------------
# Integration
#------------

@testset "ChanceConstraint integration" begin
    FactorGraph()
    @RV x ~ GaussianMeanVariance(0.0, 1.0)
    ChanceConstraint(x; G=(0.0, Inf), epsilon=0.05)
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)

    @test occursin("ruleSPChanceConstraintOutG(messages[1], (0.0, Inf), 0.05)", algo_code)
end

end # module