module Abstract_distTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, unsafeLogMean, unsafeMeanCov, unsafeMirroredLogMean, dims

@testset "Abstract_dist ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=5.2) == ProbabilityDistribution{Univariate, Abstract_dist}(Dict(:m=>0.0, :v=>1.0, :f=>5.2))
    @test_throws Exception ProbabilityDistribution(Multivariate, Abstract_dist)
end

@testset "dims" begin
    f_dummy(x) = x
    @test dims(ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=f_dummy)) == 1
end

@testset "unsafe mean and variance" begin
    f_dummy(x) = x
    sigmoid(x) = 1/(1+exp(-x))
    @test abs(unsafeMean(ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=f_dummy)) - 0.0) < 0.1
    @test abs(unsafeVar(ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=f_dummy)) - 1.0) < 0.1
    @test abs(unsafeLogMean(ProbabilityDistribution(Univariate, Abstract_dist, m=20000.0, v=0.01, f=sigmoid)) - 0.0) < 0.1
    @test abs(unsafeMirroredLogMean(ProbabilityDistribution(Univariate, Abstract_dist, m=-10000.0, v=0.1, f=sigmoid)) - 0.0) < 0.1
end

@testset "prod!" begin
    @test ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=1) * ProbabilityDistribution(Bernoulli, p=0.8) == ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=1)
    @test ProbabilityDistribution(Bernoulli, p=0.8) * ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=1) == ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=1)
end

end
