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
    @test outboundType(VBGaussianMixtureM1) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMixtureM1, [Univariate, Void, Univariate, Univariate, Univariate, Univariate]) 

    @test ruleVBGaussianMixtureM1(Univariate(Gaussian, m=8.5, v=0.5), nothing, Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=0.0, v=1.0), Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gaussian, m=8.5, v=10.0))
end

@testset "VBGaussianMixtureW1" begin
    @test VBGaussianMixtureW1 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW1) == Message{Univariate{Gamma}}
    @test isApplicable(VBGaussianMixtureW1, [Univariate, Univariate, Void, Univariate, Univariate, Univariate]) 

    @test ruleVBGaussianMixtureW1(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=5.0, v=2.0), nothing, Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gamma, a=1.1, b=1.475))
end

@testset "VBGaussianMixtureM2" begin
    @test VBGaussianMixtureM2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM2) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMixtureM2, [Univariate, Univariate, Univariate, Void, Univariate, Univariate]) 

    @test ruleVBGaussianMixtureM2(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=0.0, v=1.0), Univariate(Gamma, a=1.0, b=2.0), nothing, Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gaussian, m=8.5, v=0.625))
end

@testset "VBGaussianMixtureW2" begin
    @test VBGaussianMixtureW2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW2) == Message{Univariate{Gamma}}
    @test isApplicable(VBGaussianMixtureW2, [Univariate, Univariate, Univariate, Univariate, Void, Univariate]) 

    @test ruleVBGaussianMixtureW2(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), nothing, Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gamma, a=1.4, b=2.3000000000000003))
end

@testset "VBGaussianMixtureZ" begin
    @test VBGaussianMixtureZ <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZ) == Message{Univariate{Bernoulli}}
    @test isApplicable(VBGaussianMixtureZ, [Univariate, Univariate, Univariate, Univariate, Univariate, Void]) 

    @test ruleVBGaussianMixtureZ(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0), nothing) == Message(Univariate(Bernoulli, p=0.7713458788198754))
end

@testset "VBGaussianMixtureOut" begin
    @test VBGaussianMixtureOut <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureOut) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMixtureOut, [Void, Univariate, Univariate, Univariate, Univariate, Univariate]) 

    @test ruleVBGaussianMixtureOut(nothing, Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0), Univariate(Bernoulli, p=0.2)) == Message(Univariate(Gaussian, m=9.705882352941178, v=0.5882352941176471))
end

end #module