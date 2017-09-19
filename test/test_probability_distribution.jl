module ProbabilityDistributionTest

using Base.Test
import ForneyLab: ProbabilityDistribution, GaussianMeanVariance, PointMass, Equality, mean, var, isProper

@testset "ProbabilityDistribution" begin
    # ProbabilityDistribution should be parameterized on a node type (distribution family)
    gaussian = ProbabilityDistribution(GaussianMeanVariance, m=0.0, v=1.0)
    @test isa(gaussian, ProbabilityDistribution{GaussianMeanVariance})
    @test gaussian.params == Dict(:m=>0.0, :v=>1.0)

    # PointMass should be defined as a special family 
    point_mass = ProbabilityDistribution(PointMass, m=0.0)
    @test isa(point_mass, ProbabilityDistribution{PointMass})
    @test point_mass.params == Dict(:m=>0.0)
    @test mean(point_mass) == 0.0
    @test_throws Exception var(point_mass)
    @test isProper(point_mass)

    # Remaining DeltaFactors should not be allowed as a distribution family
    @test_throws Exception ProbabilityDistribution(Equality)
    @test_throws Exception ProbabilityDistribution(Constant, m=0.0)
end

end #module