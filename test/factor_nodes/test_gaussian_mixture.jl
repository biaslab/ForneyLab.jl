module GaussianMixtureTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: VBGaussianMixtureM1, VBGaussianMixtureW1, VBGaussianMixtureM2, VBGaussianMixtureW2, VBGaussianMixtureZ, VBGaussianMixtureOut

#-------------
# Update rules
#-------------

@testset "VBGaussianMixtureM1" begin
    @test VBGaussianMixtureM1 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM1) == Message{Gaussian}
    @test isApplicable(VBGaussianMixtureM1, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureM1(Univariate(Gaussian, m=8.5, v=0.5), nothing, Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=0.0, v=1.0), Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gaussian, m=8.5, v=10.0))
end

@testset "VBGaussianMixtureW1" begin
    @test VBGaussianMixtureW1 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW1) == Message{Gamma}
    @test isApplicable(VBGaussianMixtureW1, [ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate}, Void, ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate}]) 
    @test !isApplicable(VBGaussianMixtureW1, [ProbabilityDistribution{Multivariate}, ProbabilityDistribution{Multivariate}, Void, ProbabilityDistribution{Multivariate}, ProbabilityDistribution{MatrixVariate}, ProbabilityDistribution{Univariate}]) 

    @test ruleVBGaussianMixtureW1(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=5.0, v=2.0), nothing, Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gamma, a=1.1, b=1.475))
end

@testset "VBGaussianMixtureM2" begin
    @test VBGaussianMixtureM2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM2) == Message{Gaussian}
    @test isApplicable(VBGaussianMixtureM2, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureM2(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=0.0, v=1.0), Univariate(Gamma, a=1.0, b=2.0), nothing, Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gaussian, m=8.5, v=0.625))
end

@testset "VBGaussianMixtureW2" begin
    @test VBGaussianMixtureW2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW2) == Message{Gamma}
    @test isApplicable(VBGaussianMixtureW2, [ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate}, ProbabilityDistribution{Univariate}, Void, ProbabilityDistribution{Univariate}]) 
    @test !isApplicable(VBGaussianMixtureW2, [ProbabilityDistribution{Multivariate}, ProbabilityDistribution{Multivariate}, ProbabilityDistribution{MatrixVariate}, ProbabilityDistribution{Multivariate}, Void, ProbabilityDistribution{Univariate}]) 

    @test ruleVBGaussianMixtureW2(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), nothing, Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gamma, a=1.4, b=2.3000000000000003))
end

@testset "VBGaussianMixtureZ" begin
    @test VBGaussianMixtureZ <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZ) == Message{Bernoulli}
    @test isApplicable(VBGaussianMixtureZ, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMixtureZ(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0), nothing) == Message(Univariate(Bernoulli, p=0.7713458788198754))
end

@testset "VBGaussianMixtureOut" begin
    @test VBGaussianMixtureOut <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureOut) == Message{Gaussian}
    @test isApplicable(VBGaussianMixtureOut, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureOut(nothing, Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gaussian, m=9.705882352941178, v=0.5882352941176471))
end

end #module