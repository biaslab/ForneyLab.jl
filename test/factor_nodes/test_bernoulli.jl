module BernoulliTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPBernoulliPV, prod!, unsafeMean, unsafeVar

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Bernoulli, p=0.2)) == 0.2
    @test unsafeVar(ProbabilityDistribution(Bernoulli, p=0.5)) == 0.25
end

@testset "prod!" begin
    @test ProbabilityDistribution(Bernoulli, p=0.2) * ProbabilityDistribution(Bernoulli, p=0.8) == ProbabilityDistribution(Bernoulli, p=0.5000000000000001)
    @test_throws Exception ProbabilityDistribution(Bernoulli, p=0.0) * ProbabilityDistribution(Bernoulli, p=1.0)
end

#-------------
# Update rules
#-------------

@testset "SPBernoulliPV" begin
    @test SPBernoulliPV <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliPV) == Message{Bernoulli}
    @test isApplicable(SPBernoulliPV, [Message{PointMass}, Void]) 

    @test ruleSPBernoulliPV(Message(PointMass, m=0.2), nothing) == Message(Bernoulli, p=0.2)
end

end # module