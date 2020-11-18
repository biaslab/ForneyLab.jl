module PointMassConstraintTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPPointMassConstraintOut


#-------------
# Update rules
#-------------

@testset "SPPointMassConstraintOut" begin
    @test SPPointMassConstraintOut <: SumProductRule{PointMassConstraint}
    @test outboundType(SPPointMassConstraintOut) == Message{PointMass}
    @test isApplicable(SPPointMassConstraintOut, [Nothing]) 

    @test ruleSPPointMassConstraintOut(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == Message(Univariate, PointMass, m=0.0)
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

    @test occursin("ruleSPPointMassConstraintOut(messages[1])", algo_code)
end

end # module