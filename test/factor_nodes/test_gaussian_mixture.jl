module GaussianMixtureTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: VBGaussianMixtureM, VBGaussianMixtureW, VBGaussianMixtureZBer, VBGaussianMixtureZCat, VBGaussianMixtureOut

#-------------
# Update rules
#-------------

@testset "VBGaussianMixtureM" begin
    @test VBGaussianMixtureM <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureM) == Message{Gaussian{Precision}}
    @test !isApplicable(VBGaussianMixtureM, [Nothing, Distribution, Distribution, Distribution, Distribution, Distribution]) 
    @test !isApplicable(VBGaussianMixtureM, [Distribution, Distribution, Distribution, Nothing, Distribution, Distribution]) 
    @test !isApplicable(VBGaussianMixtureM, [Distribution, Distribution, Nothing, Distribution, Distribution]) 
    @test isApplicable(VBGaussianMixtureM, [Distribution, Distribution, Nothing, Distribution, Distribution, Distribution]) 
    @test isApplicable(VBGaussianMixtureM, [Distribution, Distribution, Distribution, Distribution, Nothing, Distribution]) 
    @test isApplicable(VBGaussianMixtureM, [Distribution, Distribution, Distribution, Distribution, Distribution, Distribution, Nothing, Distribution]) 

    @test ruleVBGaussianMixtureM(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), Distribution(Univariate, Bernoulli, p=0.2), nothing, Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian{Precision}, m=8.5, w=0.1)
    @test ruleVBGaussianMixtureM(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), Distribution(Univariate, Bernoulli, p=0.2), nothing, Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, Gaussian{Precision}, m=[8.5], w=mat(0.1))
    @test ruleVBGaussianMixtureM(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), Distribution(Univariate, Bernoulli, p=0.2), Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0), Distribution(Univariate, Gamma, a=1.0, b=2.0), nothing, Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian{Precision}, m=8.5, w=1.6)
    @test ruleVBGaussianMixtureM(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), Distribution(Univariate, Bernoulli, p=0.2), Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0)), Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), nothing, Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, Gaussian{Precision}, m=[8.5], w=mat(1.6))
    @test ruleVBGaussianMixtureM(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), Distribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian{Precision}, m=8.5, w=0.1)
    @test ruleVBGaussianMixtureM(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), Distribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, Gaussian{Precision}, m=[8.5], w=mat(0.1))
end

@testset "VBGaussianMixtureW" begin
    @test VBGaussianMixtureW <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureW) == Message{Union{Gamma, Wishart}}
    @test !isApplicable(VBGaussianMixtureW, [Nothing, Distribution, Distribution, Distribution, Distribution, Distribution]) 
    @test !isApplicable(VBGaussianMixtureW, [Distribution, Distribution, Nothing, Distribution, Distribution, Distribution]) 
    @test !isApplicable(VBGaussianMixtureW, [Distribution, Distribution, Distribution, Nothing, Distribution]) 
    @test isApplicable(VBGaussianMixtureW, [Distribution, Distribution, Distribution, Nothing, Distribution, Distribution]) 
    @test isApplicable(VBGaussianMixtureW, [Distribution, Distribution, Distribution, Distribution, Distribution, Nothing]) 
    @test isApplicable(VBGaussianMixtureW, [Distribution, Distribution, Distribution, Distribution, Distribution, Distribution, Distribution, Nothing]) 

    @test ruleVBGaussianMixtureW(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), Distribution(Univariate, Bernoulli, p=0.2), Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), nothing, Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gamma, a=1.1, b=1.475)
    @test ruleVBGaussianMixtureW(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), Distribution(Univariate, Bernoulli, p=0.2), Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), nothing, Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(MatrixVariate, Wishart, nu=2.2, v=mat(0.33898305084745767))
    @test ruleVBGaussianMixtureW(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), Distribution(Univariate, Bernoulli, p=0.2), Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), nothing) == Message(Univariate, Gamma, a=1.4, b=2.3000000000000003)
    @test ruleVBGaussianMixtureW(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), Distribution(Univariate, Bernoulli, p=0.2), Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), nothing) == Message(MatrixVariate, Wishart, nu=2.8, v=mat(0.2173913043478261))
    @test ruleVBGaussianMixtureW(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), Distribution(Univariate, Categorical, p=[0.2, 0.8]), Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), nothing, Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gamma, a=1.1, b=1.475)
    @test ruleVBGaussianMixtureW(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), Distribution(Univariate, Categorical, p=[0.2, 0.8]), Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), nothing, Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(MatrixVariate, Wishart, nu=2.2, v=mat(0.33898305084745767))
end

@testset "VBGaussianMixtureZBer" begin
    @test VBGaussianMixtureZBer <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZBer) == Message{Bernoulli}
    @test !isApplicable(VBGaussianMixtureZBer, [Distribution, Nothing, Distribution, Distribution, Distribution, Distribution, Distribution, Distribution]) 
    @test isApplicable(VBGaussianMixtureZBer, [Distribution, Nothing, Distribution, Distribution, Distribution, Distribution]) 

    @test ruleVBGaussianMixtureZBer(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), nothing, Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.7713458788198754)
    @test ruleVBGaussianMixtureZBer(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), nothing, Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Univariate, Bernoulli, p=0.7713458788198754)
end

@testset "VBGaussianMixtureZCat" begin
    @test VBGaussianMixtureZCat <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureZCat) == Message{Categorical}
    @test isApplicable(VBGaussianMixtureZCat, [Distribution, Nothing, Distribution, Distribution, Distribution, Distribution, Distribution, Distribution]) 
    @test !isApplicable(VBGaussianMixtureZCat, [Distribution, Nothing, Distribution, Distribution, Distribution, Distribution]) 

    @test ruleVBGaussianMixtureZCat(Distribution(Univariate, Gaussian{Moments}, m=8.5, v=0.5), nothing, Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Categorical, p=[0.7713458788198754, 0.22865412118012463])
    @test ruleVBGaussianMixtureZCat(Distribution(Multivariate, Gaussian{Moments}, m=[8.5], v=mat(0.5)), nothing, Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Univariate, Categorical, p=[0.7713458788198754, 0.22865412118012463])
end

@testset "VBGaussianMixtureOut" begin
    @test VBGaussianMixtureOut <: NaiveVariationalRule{GaussianMixture}
    @test outboundType(VBGaussianMixtureOut) == Message{Gaussian{Canonical}}
    @test isApplicable(VBGaussianMixtureOut, [Nothing, Distribution, Distribution, Distribution, Distribution, Distribution]) 
    @test isApplicable(VBGaussianMixtureOut, [Nothing, Distribution, Distribution, Distribution, Distribution, Distribution, Distribution, Distribution]) 

    @test ruleVBGaussianMixtureOut(nothing, Distribution(Univariate, Bernoulli, p=0.2), Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian{Canonical}, xi=16.5, w=1.7)
    @test ruleVBGaussianMixtureOut(nothing, Distribution(Univariate, Bernoulli, p=0.2), Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, Gaussian{Canonical}, xi=[16.5], w=mat(1.7))
    @test ruleVBGaussianMixtureOut(nothing, Distribution(Univariate, Categorical, p=[0.2, 0.8]), Distribution(Univariate, Gaussian{Moments}, m=5.0, v=2.0), Distribution(Univariate, Gamma, a=1.0, b=2.0), Distribution(Univariate, Gaussian{Moments}, m=10.0, v=3.0), Distribution(Univariate, Gamma, a=2.0, b=1.0)) == Message(Univariate, Gaussian{Canonical}, xi=16.5, w=1.7)
    @test ruleVBGaussianMixtureOut(nothing, Distribution(Univariate, Categorical, p=[0.2, 0.8]), Distribution(Multivariate, Gaussian{Moments}, m=[5.0], v=mat(2.0)), Distribution(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Distribution(Multivariate, Gaussian{Moments}, m=[10.0], v=mat(3.0)), Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(0.5))) == Message(Multivariate, Gaussian{Canonical}, xi=[16.5], w=mat(1.7))
end

@testset "averageEnergy" begin
    # Univariate
    marg_out = Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    marg_switch = Distribution(Univariate, Bernoulli, p=0.2)
    marg_mean_1 = Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0)
    marg_prec_1 = Distribution(Univariate, Gamma, a=2.0, b=3.0)
    marg_mean_2 = Distribution(Univariate, Gaussian{Moments}, m=3.0, v=4.0)
    marg_prec_2 = Distribution(Univariate, Gamma, a=4.0, b=5.0)
    ref_val = 0.2*averageEnergy(Gaussian{Precision}, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(Gaussian{Precision}, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    marg_switch = Distribution(Univariate, Categorical, p=[0.2, 0.8])
    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    # Multivariate
    marg_out = Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0))
    marg_switch = Distribution(Univariate, Bernoulli, p=0.2)
    marg_mean_1 = Distribution(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0))
    marg_prec_1 = Distribution(MatrixVariate, Wishart, nu=4.0, v=mat(1/6))
    marg_mean_2 = Distribution(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(4.0))
    marg_prec_2 = Distribution(MatrixVariate, Wishart, nu=8.0, v=mat(0.1))
    ref_val = 0.2*averageEnergy(Gaussian{Precision}, marg_out, marg_mean_1, marg_prec_1) +
              0.8*averageEnergy(Gaussian{Precision}, marg_out, marg_mean_2, marg_prec_2)

    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val

    marg_switch = Distribution(Univariate, Categorical, p=[0.2, 0.8])
    @test averageEnergy(GaussianMixture, marg_out, marg_switch, marg_mean_1, marg_prec_1, marg_mean_2, marg_prec_2) == ref_val
end

end #module