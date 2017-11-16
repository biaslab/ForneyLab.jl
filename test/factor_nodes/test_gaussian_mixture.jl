module GaussianMixtureTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: VBGaussianMixtureM, VBGaussianMixtureW, VBGaussianMixtureZBer, VBGaussianMixtureZCat, VBGaussianMixtureOut

#-------------
# Update rules
#-------------

@testset "VBGaussianMixtureM" begin
    @test VBGaussianMixtureM <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM) == Message{Gaussian}
    @test !isApplicable(VBGaussianMixtureM, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), nothing, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian, m=8.5, v=10.0)
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), nothing, ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate, Gaussian, m=[8.5], v=[9.999999999999998].')
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), nothing, ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian, m=8.5, v=0.625)
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), nothing, ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate, Gaussian, m=[8.5], v=[0.6249999999999999].')
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian, m=8.5, v=10.0)
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate, Gaussian, m=[8.5], v=[9.999999999999998].')
end

@testset "VBGaussianMixtureW" begin
    @test VBGaussianMixtureW <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW) == Message{Scale}
    @test !isApplicable(VBGaussianMixtureW, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void]) 
    @test isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), nothing, ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gamma, a=1.1, b=1.475)
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), nothing, ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(MatrixVariate, Wishart, nu=2.2, v=[0.33898305084745767].')
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), nothing) == Message(Univariate, Gamma, a=1.4, b=2.3000000000000003)
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), nothing) == Message(MatrixVariate, Wishart, nu=2.8, v=[0.2173913043478261].')
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), nothing, ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gamma, a=1.1, b=1.475)
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), nothing, ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(MatrixVariate, Wishart, nu=2.2, v=[0.33898305084745767].')
end

@testset "VBGaussianMixtureZBer" begin
    @test VBGaussianMixtureZBer <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZBer) == Message{Bernoulli}
    @test !isApplicable(VBGaussianMixtureZBer, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureZBer, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureZBer(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), nothing, ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.7713458788198754)
    @test ruleVBGaussianMixtureZBer(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), nothing, ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Univariate, Bernoulli, p=0.7713458788198754)
end

@testset "VBGaussianMixtureZCat" begin
    @test VBGaussianMixtureZCat <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZCat) == Message{Categorical}
    @test isApplicable(VBGaussianMixtureZCat, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureZCat, [ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureZCat(ProbabilityDistribution(Univariate, Gaussian, m=8.5, v=0.5), nothing, ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Categorical, p=[0.7713458788198754, 0.2286541211801246])
    @test ruleVBGaussianMixtureZCat(ProbabilityDistribution(Multivariate, Gaussian, m=[8.5], v=[0.5].'), nothing, ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Univariate, Categorical, p=[0.7713458788198754, 0.2286541211801246])
end

@testset "VBGaussianMixtureOut" begin
    @test VBGaussianMixtureOut <: VariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureOut) == Message{Gaussian}
    @test isApplicable(VBGaussianMixtureOut, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureOut, [Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian, m=9.705882352941178, v=0.5882352941176471)
    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate, Gaussian, m=[9.705882352941176], v=[0.588235294117647].')
    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, Gaussian, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, Gaussian, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian, m=9.705882352941178, v=0.5882352941176471)
    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Multivariate, Gaussian, m=[5.0], v=[2.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), ProbabilityDistribution(Multivariate, Gaussian, m=[10.0], v=[3.0].'), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[0.5].')) == Message(Multivariate, Gaussian, m=[9.705882352941176], v=[0.588235294117647].')
end

@testset "averageEnergy" begin
    # Univariate
    marg_out = ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)
    marg_switch = ProbabilityDistribution(Univariate, Bernoulli, p=0.2)
    marg_mean_1 = ProbabilityDistribution(Univariate, Gaussian, m=1.0, v=2.0)
    marg_prec_1 = ProbabilityDistribution(Univariate, Gamma, a=2.0, b=3.0)
    marg_mean_2 = ProbabilityDistribution(Univariate, Gaussian, m=3.0, v=4.0)
    marg_prec_2 = ProbabilityDistribution(Univariate, Gamma, a=4.0, b=5.0)
    ref_val = 0.2*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    marg_switch = ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8])
    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    # Multivariate
    marg_out = ProbabilityDistribution(Multivariate, Gaussian, m=[0.0], v=[1.0].')
    marg_switch = ProbabilityDistribution(Univariate, Bernoulli, p=0.2)
    marg_mean_1 = ProbabilityDistribution(Multivariate, Gaussian, m=[1.0], v=[2.0].')
    marg_prec_1 = ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=[1/6].')
    marg_mean_2 = ProbabilityDistribution(Multivariate, Gaussian, m=[3.0], v=[4.0].')
    marg_prec_2 = ProbabilityDistribution(MatrixVariate, Wishart, nu=8.0, v=[0.1].')
    ref_val = 0.2*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    marg_switch = ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8])
    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val
end

end #module