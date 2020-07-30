module PointMassConstraintTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPPointMassConstraintOutG


#-------------
# Update rules
#-------------

@testset "SPPointMassConstraintOutG" begin
    @test SPPointMassConstraintOutG <: SumProductRule{PointMassConstraint}
    @test outboundType(SPPointMassConstraintOutG) == Message{PointMass}
    @test isApplicable(SPPointMassConstraintOutG, [Nothing]) 

    @test ruleSPPointMassConstraintOutG(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == Message(Univariate, PointMass, m=0.0)
end


#------------
# Integration
#------------

@testset "PointMassConstraint integration" begin
    fg = FactorGraph()
    @RV x ~ GaussianMeanVariance(0.0, 1.0)
    PointMassConstraint(x)
    
    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    algo_code = algorithmSourceCode(algo)

    @test occursin("ruleSPPointMassConstraintOutG(messages[1])", algo_code)
end

end # module