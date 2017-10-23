module GaussianMixtureTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, VBGaussianMixture1, VBGaussianMixture2, VBGaussianMixture3, VBGaussianMixture4, VBGaussianMixture5, VBGaussianMixture6


#-------------
# Update rules
#-------------

@testset "VBGaussianMixture1" begin
    @test VBGaussianMixture1 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixture1) == Message{Gaussian}
    @test isApplicable(VBGaussianMixture1, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixture1(nothing, ProbabilityDistribution(Gamma, a=1.0, b=2.0), ProbabilityDistribution(Gaussian, m=0.0, v=1.0), ProbabilityDistribution(Gamma, a=2.0, b=1.0), ProbabilityDistribution(Bernoulli, p=0.2), ProbabilityDistribution(Gaussian, m=8.5, v=0.5)) == Message(Gaussian, m=8.5, v=10.0)
end

@testset "VBGaussianMixture2" begin
    @test VBGaussianMixture2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixture2) == Message{Gamma}
    @test isApplicable(VBGaussianMixture2, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixture2(ProbabilityDistribution(Gaussian, m=5.0, v=2.0), nothing, ProbabilityDistribution(Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Gamma, a=2.0, b=1.0), ProbabilityDistribution(Bernoulli, p=0.2), ProbabilityDistribution(Gaussian, m=8.5, v=0.5)) == Message(Gamma, a=1.1, b=1.475)
end

@testset "VBGaussianMixture3" begin
    @test VBGaussianMixture3 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixture3) == Message{Gaussian}
    @test isApplicable(VBGaussianMixture3, [ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixture3(ProbabilityDistribution(Gaussian, m=0.0, v=1.0), ProbabilityDistribution(Gamma, a=1.0, b=2.0), nothing, ProbabilityDistribution(Gamma, a=2.0, b=1.0), ProbabilityDistribution(Bernoulli, p=0.2), ProbabilityDistribution(Gaussian, m=8.5, v=0.5)) == Message(Gaussian, m=8.5, v=0.625)
end

@testset "VBGaussianMixture4" begin
    @test VBGaussianMixture4 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixture4) == Message{Gamma}
    @test isApplicable(VBGaussianMixture4, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixture4(ProbabilityDistribution(Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Gamma, a=1.0, b=2.0), ProbabilityDistribution(Gaussian, m=10.0, v=3.0), nothing, ProbabilityDistribution(Bernoulli, p=0.2), ProbabilityDistribution(Gaussian, m=8.5, v=0.5)) == Message(Gamma, a=1.4, b=2.3000000000000003)
end

@testset "VBGaussianMixture5" begin
    @test VBGaussianMixture5 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixture5) == Message{Bernoulli}
    @test isApplicable(VBGaussianMixture5, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution]) 

    @test ruleVBGaussianMixture5(ProbabilityDistribution(Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Gamma, a=1.0, b=2.0), ProbabilityDistribution(Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Gamma, a=2.0, b=1.0), nothing, ProbabilityDistribution(Gaussian, m=8.5, v=0.5)) == Message(Bernoulli, p=0.7713458788198754)
end

@testset "VBGaussianMixture6" begin
    @test VBGaussianMixture6 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixture6) == Message{Gaussian}
    @test isApplicable(VBGaussianMixture6, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMixture6(ProbabilityDistribution(Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Gamma, a=1.0, b=2.0), ProbabilityDistribution(Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Gamma, a=2.0, b=1.0), ProbabilityDistribution(Bernoulli, p=0.2), nothing) == Message(Gaussian, m=9.705882352941178, v=0.5882352941176471)
end

end #module