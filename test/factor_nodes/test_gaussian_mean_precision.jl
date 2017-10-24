module GaussianMeanPrecisionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, VBGaussianMeanPrecision1, VBGaussianMeanPrecision2, VBGaussianMeanPrecision3


#-------------
# Update rules
#-------------

@testset "VBGaussianMeanPrecision1" begin
    @test VBGaussianMeanPrecision1 <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecision1) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanPrecision1, [Void, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMeanPrecision1, [ProbabilityDistribution, Void, ProbabilityDistribution]) 

    @test ruleVBGaussianMeanPrecision1(nothing, ProbabilityDistribution(Gamma, a=1.0, b=2.0), ProbabilityDistribution(Gaussian, m=3.0, v=4.0)) == Message(Gaussian, m=3.0, w=0.5)
end

@testset "VBGaussianMeanPrecision2" begin
    @test VBGaussianMeanPrecision2 <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecision2) == Message{Gamma}
    @test isApplicable(VBGaussianMeanPrecision2, [ProbabilityDistribution, Void, ProbabilityDistribution]) 

    @test ruleVBGaussianMeanPrecision2(ProbabilityDistribution(Gaussian, m=1.0, v=2.0), Void, ProbabilityDistribution(Gaussian, m=3.0, v=4.0)) == Message(Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2))
end

@testset "VBGaussianMeanPrecision3" begin
    @test VBGaussianMeanPrecision3 <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecision3) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanPrecision3, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMeanPrecision3(ProbabilityDistribution(Gaussian, m=3.0, v=4.0), ProbabilityDistribution(Gamma, a=1.0, b=2.0), nothing) == Message(Gaussian, m=3.0, w=0.5)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Gaussian, m=0.0, w=2.0)) == averageEnergy(GaussianMeanPrecision, Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0), Univariate(Gaussian, m=0.0, w=2.0))
    @test differentialEntropy(Univariate(Gaussian, m=0.0, w=2.0)) == differentialEntropy(Multivariate(Gaussian, m=[0.0], w=[2.0].'))
    @test averageEnergy(GaussianMeanPrecision, Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0), Univariate(Gaussian, m=0.0, w=2.0)) == averageEnergy(GaussianMeanPrecision, Multivariate(PointMass, m=[0.0]), Multivariate(PointMass, m=[2.0]), Multivariate(Gaussian, m=[0.0], W=[2.0].'))
end

end #module