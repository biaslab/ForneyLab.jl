module GaussianMeanPrecisionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: VBGaussianMeanPrecisionOut, VBGaussianMeanPrecisionM, VBGaussianMeanPrecisionW


#-------------
# Update rules
#-------------

@testset "VBGaussianMeanPrecisionM" begin
    @test VBGaussianMeanPrecisionM <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionM) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanPrecisionM, [ProbabilityDistribution, Void, ProbabilityDistribution]) 
    @test !isApplicable(VBGaussianMeanPrecisionM, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMeanPrecisionM(Univariate(Gaussian, m=3.0, v=4.0), nothing, Univariate(Gamma, a=1.0, b=2.0)) == Message(Univariate(Gaussian, m=3.0, w=0.5))
    @test ruleVBGaussianMeanPrecisionM(Multivariate(Gaussian, m=[3.0], v=[4.0].'), nothing, MatrixVariate(Wishart, v=[0.25].', nu=2.0)) == Message(Multivariate(Gaussian, m=[3.0], w=[0.5].'))
end

@testset "VBGaussianMeanPrecisionW" begin
    @test VBGaussianMeanPrecisionW <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionW) == Message{Scale}
    @test isApplicable(VBGaussianMeanPrecisionW, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBGaussianMeanPrecisionW(Univariate(Gaussian, m=3.0, v=4.0), Univariate(Gaussian, m=1.0, v=2.0), nothing) == Message(Univariate(Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2)))
    @test ruleVBGaussianMeanPrecisionW(Multivariate(Gaussian, m=[3.0], v=[4.0].'), Multivariate(Gaussian, m=[1.0], v=[2.0].'), nothing) == Message(MatrixVariate(Wishart, v=[1.0/(2.0 + 4.0 + (3.0 - 1.0)^2)].', nu=3.0))
end

@testset "VBGaussianMeanPrecisionOut" begin
    @test VBGaussianMeanPrecisionOut <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionOut) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanPrecisionOut, [Void, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBGaussianMeanPrecisionOut(nothing, Univariate(Gaussian, m=3.0, v=4.0), Univariate(Gamma, a=1.0, b=2.0)) == Message(Univariate(Gaussian, m=3.0, w=0.5))
    @test ruleVBGaussianMeanPrecisionOut(nothing, Multivariate(Gaussian, m=[3.0], v=[4.0].'), MatrixVariate(Wishart, v=[0.25].', nu=2.0)) == Message(Multivariate(Gaussian, m=[3.0], w=[0.5].'))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Gaussian, m=0.0, w=2.0)) == averageEnergy(GaussianMeanPrecision, Univariate(Gaussian, m=0.0, w=2.0), Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0))
    @test differentialEntropy(Univariate(Gaussian, m=0.0, w=2.0)) == differentialEntropy(Multivariate(Gaussian, m=[0.0], w=[2.0].'))
    @test averageEnergy(GaussianMeanPrecision, Univariate(Gaussian, m=0.0, w=2.0), Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0)) == averageEnergy(GaussianMeanPrecision, Multivariate(Gaussian, m=[0.0], w=[2.0].'), Multivariate(PointMass, m=[0.0]), MatrixVariate(PointMass, m=[2.0].'))
end

end #module