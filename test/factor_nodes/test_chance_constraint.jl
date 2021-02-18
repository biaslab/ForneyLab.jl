module ChanceConstraintTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPChanceConstraintOutG, truncatedGaussianMoments


#--------
# Helpers
#--------

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

    @test ruleSPChanceConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), (0.0, Inf), 0.05) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.9470311004635064, w=2.2101933871827293)
    @test ruleSPChanceConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), (-Inf, Inf), 0.05) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
end


#------------
# Integration
#------------

@testset "ChanceConstraint integration" begin
    fg = FactorGraph()
    @RV x ~ GaussianMeanVariance(0.0, 1.0)
    ChanceConstraint(x; G=(0.0, Inf), epsilon=0.05)
    
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    algo_code = algorithmSourceCode(algo)

    @test occursin("ruleSPChanceConstraintOutG(messages[1], (0.0, Inf), 0.05)", algo_code)
end

end # module