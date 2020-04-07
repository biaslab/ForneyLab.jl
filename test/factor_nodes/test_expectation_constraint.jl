module ExpectationConstraintTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPExpectationConstraintOutG
using ForneyLab: expectedValue, constrainedMoments

f(x::Float64) = x # Trivial
g(x::Float64) = 1.0*(x > 0.0) # Non-trivial


#--------
# Helpers
#--------

@testset "expectedValue" begin
    @test expectedValue(2.0, 0.0, 1.0, f) == 1.999999026302908
    @test expectedValue(2.0, 0.0, 1.0, g) == 0.8807970779778822
end

@testset "constrainedMoments" begin
    @test constrainedMoments(2.0, 0.0, 1.0, f) == (1.999999026302908, 0.9999990688765581)
    @test constrainedMoments(2.0, 0.0, 1.0, g) == (0.6155633430596962, 0.6210829316309023)
end


#-------------
# Update rules
#-------------

@testset "SPExpectationConstraintOutG" begin
    @test SPExpectationConstraintOutG <: SumProductRule{ExpectationConstraint}
    @test outboundType(SPExpectationConstraintOutG) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPExpectationConstraintOutG, [Nothing]) 

    @test ruleSPExpectationConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), g, 0.95, 0.5) == Message(Univariate, GaussianWeightedMeanPrecision, xi=1.5449452063991789, w=1.123838185932649)
    @test ruleSPExpectationConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), f, 2.0, 0.0) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0000018622402607, w=9.31120130553964e-7)
end


#------------
# Integration
#------------

@testset "ExpectationConstraint integration" begin
    FactorGraph()
    @RV x ~ GaussianMeanVariance(0.0, 1.0)
    ExpectationConstraint(x; g=f, G=1.0, eta_init=0.0)
    algo = sumProductAlgorithm(x)
    algo_code = algorithmSourceCode(algo)

    @test occursin("ruleSPExpectationConstraintOutG(messages[1], f, 1.0, 0.0)", algo_code)
end

end # module