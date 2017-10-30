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
    @test outboundType(VBGaussianMeanPrecisionM) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMeanPrecisionM, [Univariate, Void, Univariate]) 
    @test !isApplicable(VBGaussianMeanPrecisionM, [Univariate, Univariate, Void]) 

    @test ruleVBGaussianMeanPrecisionM(Univariate(Gaussian, m=3.0, v=4.0), nothing, Univariate(Gamma, a=1.0, b=2.0)) == Message(Univariate(Gaussian, m=3.0, w=0.5))
end

@testset "VBGaussianMeanPrecisionW" begin
    @test VBGaussianMeanPrecisionW <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionW) == Message{Univariate{Gamma}}
    @test isApplicable(VBGaussianMeanPrecisionW, [Univariate, Univariate, Void]) 

    @test ruleVBGaussianMeanPrecisionW(Univariate(Gaussian, m=3.0, v=4.0), Univariate(Gaussian, m=1.0, v=2.0), nothing) == Message(Univariate(Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2)))
end

@testset "VBGaussianMeanPrecisionOut" begin
    @test VBGaussianMeanPrecisionOut <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionOut) == Message{Univariate{Gaussian}}
    @test isApplicable(VBGaussianMeanPrecisionOut, [Void, Univariate, Univariate]) 

    @test ruleVBGaussianMeanPrecisionOut(nothing, Univariate(Gaussian, m=3.0, v=4.0), Univariate(Gamma, a=1.0, b=2.0)) == Message(Univariate(Gaussian, m=3.0, w=0.5))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Univariate(Gaussian, m=0.0, w=2.0)) == averageEnergy(GaussianMeanPrecision, Univariate(Gaussian, m=0.0, w=2.0), Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0))
    @test differentialEntropy(Univariate(Gaussian, m=0.0, w=2.0)) == differentialEntropy(Multivariate(Gaussian, m=[0.0], w=[2.0].'))
    @test averageEnergy(GaussianMeanPrecision, Univariate(Gaussian, m=0.0, w=2.0), Univariate(PointMass, m=0.0), Univariate(PointMass, m=2.0)) == averageEnergy(GaussianMeanPrecision, Multivariate(Gaussian, m=[0.0], w=[2.0].'), Multivariate(PointMass, m=[0.0]), MatrixVariate(PointMass, m=[2.0].'))
end

end #module