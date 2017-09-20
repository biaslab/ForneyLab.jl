module GaussianWeightedMeanPrecisionTest

using Base.Test
using ForneyLab
import ForneyLab: ensureParameters!, prod!, ==

@testset "ensureParameters!" begin
    dist = ProbabilityDistribution(Gaussian, m=0.0, v=1.0)
    @test ensureParameters!(dist, (:m, :v)).params == Dict(:m => 0.0, :v => 1.0)
    dist = ProbabilityDistribution(Gaussian, m=0.0, v=1.0)
    @test ensureParameters!(dist, (:xi, :w)).params == Dict(:m => 0.0, :v => 1.0, :xi => 0.0, :w => 1.0)
end

@testset "==" begin
    @test ProbabilityDistribution(Gaussian, xi=0.0, w=1.0) == ProbabilityDistribution(Gaussian, xi=0.0, w=1.0)
    @test ProbabilityDistribution(Gaussian, xi=0.0, w=1.0) == ProbabilityDistribution(Gaussian, m=0.0, v=1.0)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Gaussian, xi=0.0, w=1.0) * ProbabilityDistribution(Gaussian, xi=0.0, w=1.0) == ProbabilityDistribution(Gaussian, xi=0.0, w=2.0)
    @test ProbabilityDistribution(Gaussian, m=0.0, v=1.0) * ProbabilityDistribution(Gaussian, m=0.0, v=1.0) == ProbabilityDistribution(Gaussian, xi=0.0, w=2.0)
end

end #module