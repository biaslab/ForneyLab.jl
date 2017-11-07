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
    @test isApplicable(VBGaussianMixtureM1, [ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureM1(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Bernoulli, p=0.2), nothing, Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=0.0, v=1.0), Univariate(Gamma, a=2.0, b=1.0)) == Message(Univariate(Gaussian, m=8.5, v=10.0))
    @test ruleVBGaussianMixtureM1(Multivariate(Gaussian, m=[8.5], v=[0.5].'), Univariate(Bernoulli, p=0.2), nothing, MatrixVariate(Wishart, nu=2.0, v=[0.25].'), Multivariate(Gaussian, m=[0.0], v=[1.0].'), MatrixVariate(Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate(Gaussian, m=[8.5], v=[9.999999999999998].'))
end

@testset "VBGaussianMixtureW1" begin
    @test VBGaussianMixtureW1 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW1) == Message{AbstractGamma}
    @test isApplicable(VBGaussianMixtureW1, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureW1(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Bernoulli, p=0.2), Univariate(Gaussian, m=5.0, v=2.0), nothing, Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0)) == Message(Univariate(Gamma, a=1.1, b=1.475))
    @test ruleVBGaussianMixtureW1(Multivariate(Gaussian, m=[8.5], v=[0.5].'), Univariate(Bernoulli, p=0.2), Multivariate(Gaussian, m=[5.0], v=[2.0].'), nothing, Multivariate(Gaussian, m=[10.0], v=[3.0].'), MatrixVariate(Wishart, nu=4.0, v=[0.5].')) == Message(MatrixVariate(Wishart, nu=2.2, v=[0.33898305084745767].'))
end

@testset "VBGaussianMixtureM2" begin
    @test VBGaussianMixtureM2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM2) == Message{Gaussian}
    @test isApplicable(VBGaussianMixtureM2, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureM2(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Bernoulli, p=0.2), Univariate(Gaussian, m=0.0, v=1.0), Univariate(Gamma, a=1.0, b=2.0), nothing, Univariate(Gamma, a=2.0, b=1.0)) == Message(Univariate(Gaussian, m=8.5, v=0.625))
    @test ruleVBGaussianMixtureM2(Multivariate(Gaussian, m=[8.5], v=[0.5].'), Univariate(Bernoulli, p=0.2), Multivariate(Gaussian, m=[0.0], v=[1.0].'), MatrixVariate(Wishart, nu=2.0, v=[0.25].'), nothing, MatrixVariate(Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate(Gaussian, m=[8.5], v=[0.6249999999999999].'))
end

@testset "VBGaussianMixtureW2" begin
    @test VBGaussianMixtureW2 <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW2) == Message{AbstractGamma}
    @test isApplicable(VBGaussianMixtureW2, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMixtureW2(Univariate(Gaussian, m=8.5, v=0.5), Univariate(Bernoulli, p=0.2), Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), nothing) == Message(Univariate(Gamma, a=1.4, b=2.3000000000000003))
    @test ruleVBGaussianMixtureW2(Multivariate(Gaussian, m=[8.5], v=[0.5].'), Univariate(Bernoulli, p=0.2), Multivariate(Gaussian, m=[5.0], v=[2.0].'), MatrixVariate(Wishart, nu=2.0, v=[0.25].'), Multivariate(Gaussian, m=[10.0], v=[3.0].'), nothing) == Message(MatrixVariate(Wishart, nu=2.8, v=[0.2173913043478261].'))
end

@testset "VBGaussianMixtureZ" begin
    @test VBGaussianMixtureZ <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZ) == Message{Bernoulli}
    @test isApplicable(VBGaussianMixtureZ, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureZ(Univariate(Gaussian, m=8.5, v=0.5), nothing, Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0)) == Message(Univariate(Bernoulli, p=0.7713458788198754))
    @test ruleVBGaussianMixtureZ(Multivariate(Gaussian, m=[8.5], v=[0.5].'), nothing, Multivariate(Gaussian, m=[5.0], v=[2.0].'), MatrixVariate(Wishart, nu=2.0, v=[0.25].'), Multivariate(Gaussian, m=[10.0], v=[3.0].'), MatrixVariate(Wishart, nu=4.0, v=[0.5].')) == Message(Univariate(Bernoulli, p=0.7713458788198754))
end

@testset "VBGaussianMixtureOut" begin
    @test VBGaussianMixtureOut <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureOut) == Message{Gaussian}
    @test isApplicable(VBGaussianMixtureOut, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureOut(nothing, Univariate(Bernoulli, p=0.2), Univariate(Gaussian, m=5.0, v=2.0), Univariate(Gamma, a=1.0, b=2.0), Univariate(Gaussian, m=10.0, v=3.0), Univariate(Gamma, a=2.0, b=1.0)) == Message(Univariate(Gaussian, m=9.705882352941178, v=0.5882352941176471))
    @test ruleVBGaussianMixtureOut(nothing, Univariate(Bernoulli, p=0.2), Multivariate(Gaussian, m=[5.0], v=[2.0].'), MatrixVariate(Wishart, nu=2.0, v=[0.25].'), Multivariate(Gaussian, m=[10.0], v=[3.0].'), MatrixVariate(Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate(Gaussian, m=[9.705882352941176], v=[0.588235294117647].'))
end

@testset "averageEnergy" begin
    # Univariate
    marg_out = Univariate(Gaussian, m=0.0, v=1.0)
    marg_switch = Univariate(Bernoulli, p=0.2)
    marg_mean_1 = Univariate(Gaussian, m=1.0, v=2.0)
    marg_prec_1 = Univariate(Gamma, a=2.0, b=3.0)
    marg_mean_2 = Univariate(Gaussian, m=3.0, v=4.0)
    marg_prec_2 = Univariate(Gamma, a=4.0, b=5.0)
    ref_val = 0.2*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    # Multivariate
    marg_out = Multivariate(Gaussian, m=[0.0], v=[1.0].')
    marg_switch = Univariate(Bernoulli, p=0.2)
    marg_mean_1 = Multivariate(Gaussian, m=[1.0], v=[2.0].')
    marg_prec_1 = MatrixVariate(Wishart, nu=4.0, v=[1/6].')
    marg_mean_2 = Multivariate(Gaussian, m=[3.0], v=[4.0].')
    marg_prec_2 = MatrixVariate(Wishart, nu=8.0, v=[0.1].')
    ref_val = 0.2*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val
end

end #module