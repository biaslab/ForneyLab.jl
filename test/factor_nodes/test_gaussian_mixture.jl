module GaussianMixtureTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: VBGaussianMixtureM, VBGaussianMixtureW, VBGaussianMixtureZBer, VBGaussianMixtureZCat, VBGaussianMixtureOut

#-------------
# Update rules
#-------------

@testset "VBGaussianMixtureM" begin
    @test VBGaussianMixtureM <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM) == Message{GaussianMeanPrecision}
    @test !isApplicable(VBGaussianMixtureM, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureM, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), nothing, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, GaussianMeanPrecision, m=8.5, w=0.1)
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), nothing, ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, GaussianMeanPrecision, m=[8.5], w=mat(0.1))
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), nothing, ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, GaussianMeanPrecision, m=8.5, w=1.6)
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), nothing, ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, GaussianMeanPrecision, m=[8.5], w=mat(1.6))
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, GaussianMeanPrecision, m=8.5, w=0.1)
    @test ruleVBGaussianMixtureM(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, GaussianMeanPrecision, m=[8.5], w=mat(0.1))
end

@testset "VBGaussianMixtureW" begin
    @test VBGaussianMixtureW <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW) == Message{Union{Gamma, Wishart}}
    @test !isApplicable(VBGaussianMixtureW, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing]) 
    @test isApplicable(VBGaussianMixtureW, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing]) 

    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gamma, a=1.1, b=1.475)
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(MatrixVariate, Wishart, nu=2.2, v=mat(0.33898305084745767))
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), nothing) == Message(Univariate, Gamma, a=1.4, b=2.3000000000000003)
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), nothing) == Message(MatrixVariate, Wishart, nu=2.8, v=mat(0.2173913043478261))
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gamma, a=1.1, b=1.475)
    @test ruleVBGaussianMixtureW(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(MatrixVariate, Wishart, nu=2.2, v=mat(0.33898305084745767))
end

@testset "VBGaussianMixtureZBer" begin
    @test VBGaussianMixtureZBer <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZBer) == Message{Bernoulli}
    @test !isApplicable(VBGaussianMixtureZBer, [ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureZBer, [ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureZBer(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.7713458788198754)
    @test ruleVBGaussianMixtureZBer(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Univariate, Bernoulli, p=0.7713458788198754)
end

@testset "VBGaussianMixtureZCat" begin
    @test VBGaussianMixtureZCat <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZCat) == Message{Categorical}
    @test isApplicable(VBGaussianMixtureZCat, [ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMixtureZCat, [ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureZCat(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=8.5, v=0.5), nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Categorical, p=[0.7713458788198754, 0.22865412118012463])
    @test ruleVBGaussianMixtureZCat(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[8.5], v=mat(0.5)), nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Univariate, Categorical, p=[0.7713458788198754, 0.22865412118012463])
end

@testset "VBGaussianMixtureOut" begin
    @test VBGaussianMixtureOut <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureOut) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(VBGaussianMixtureOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 
    @test isApplicable(VBGaussianMixtureOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=16.5, w=1.7)
    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[16.5], w=mat(1.7))
    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=5.0, v=2.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=10.0, v=3.0), ProbabilityDistribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=16.5, w=1.7)
    @test ruleVBGaussianMixtureOut(nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[5.0], v=mat(2.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[10.0], v=mat(3.0)), ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[16.5], w=mat(1.7))
end

@testset "averageEnergy" begin
    # Univariate
    marg_out = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    marg_switch = ProbabilityDistribution(Univariate, Bernoulli, p=0.2)
    marg_mean_1 = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    marg_prec_1 = ProbabilityDistribution(Univariate, Gamma, a=2.0, b=3.0)
    marg_mean_2 = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0)
    marg_prec_2 = ProbabilityDistribution(Univariate, Gamma, a=4.0, b=5.0)
    ref_val = 0.2*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    marg_switch = ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8])
    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    # Multivariate
    marg_out = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    marg_switch = ProbabilityDistribution(Univariate, Bernoulli, p=0.2)
    marg_mean_1 = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
    marg_prec_1 = ProbabilityDistribution(MatrixVariate, Wishart, nu=4.0, v=mat(1/6))
    marg_mean_2 = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0))
    marg_prec_2 = ProbabilityDistribution(MatrixVariate, Wishart, nu=8.0, v=mat(0.1))
    ref_val = 0.2*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(GaussianMeanPrecision, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    marg_switch = ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8])
    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val
end

end #module