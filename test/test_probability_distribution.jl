module ProbabilityDistributionTest

using Base.Test
import ForneyLab: ProbabilityDistribution, GaussianMeanVariance, PointMass, Equality, mean, var, isProper, gaussianQuadrature

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
    @test var(point_mass) == 0.0
    @test isProper(point_mass)

    # Remaining DeltaFactors should not be allowed as a distribution family
    @test_throws Exception ProbabilityDistribution(Equality)
    @test_throws Exception ProbabilityDistribution(Clamp, m=0.0)

    @test ProbabilityDistribution(PointMass, m=0.0) == ProbabilityDistribution(PointMass, m=0.0)
    @test ProbabilityDistribution(PointMass, m=0.0) != ProbabilityDistribution(GaussianMeanVariance, m=0.0, v=1.0)
    @test ProbabilityDistribution(PointMass, m=0.0) != ProbabilityDistribution(PointMass, m=1.0)
end

@testset "gaussianQuadrature" begin
    @test gaussianQuadrature(x -> 1.0, m=-1.0, v=2.0) == 0.9999996733487053
    @test gaussianQuadrature(x -> x, m=-1.0, v=2.0) == -0.9999996733487055
    @test gaussianQuadrature(x -> x^2, m=-1.0, v=2.0) == 3.0000013419448215
end

end #module