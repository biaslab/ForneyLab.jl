module GaussianMeanPrecisionTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianMeanPrecisionOutNPP, SPGaussianMeanPrecisionMPNP, SPGaussianMeanPrecisionOutNGP, SPGaussianMeanPrecisionMGNP, VBGaussianMeanPrecisionOut, VBGaussianMeanPrecisionM, VBGaussianMeanPrecisionW, SVBGaussianMeanPrecisionOutVGD, SVBGaussianMeanPrecisionMGVD, SVBGaussianMeanPrecisionW, MGaussianMeanPrecisionGGD

@testset "dims" begin
    @test dims(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0)) == 1
    @test dims(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=ones(2), w=diageye(2))) == 2
end

@testset "vague" begin
    @test vague(GaussianMeanPrecision) == ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=tiny)
    @test vague(GaussianMeanPrecision, 2) == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=zeros(2), w=tiny*eye(2))
    @test vague(GaussianMeanPrecision, (2,)) == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=zeros(2), w=tiny*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0))
    @test !isProper(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=-1.0))

    # Multivariate
    @test isProper(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0)))
    @test isProper(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=ones(2), w=diageye(2)))
    @test !isProper(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0) == ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0)
    @test ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0) == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)

    # Multivariate
    @test ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0)) == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0))
    @test ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0)) == ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 2.0
    @test unsafeMode(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 2.0
    @test unsafeVar(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 0.25
    @test unsafeCov(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 0.25
    @test unsafeMeanCov(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == (2.0, 0.25)
    @test unsafePrecision(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 4.0
    @test unsafeWeightedMean(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 8.0
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == (8.0, 4.0)

    # Multivariate
    @test unsafeMean(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeMode(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeVar(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [0.25]
    @test unsafeCov(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == mat(0.25)
    @test unsafeMeanCov(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == ([2.0], mat(0.25))
    @test unsafePrecision(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == mat(4.0)
    @test unsafeWeightedMean(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [8.0]
    @test unsafeWeightedMeanPrecision(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == ([8.0], mat(4.0))
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=1.0, w=2.0), 1.0), -0.5723649429247)
    @test isapprox(logPdf(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[1.0, 1.0], w=[2.0 0.0; 0.0 2.0]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(ProbabilityDistribution{Univariate, GaussianMeanPrecision}, ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=8.0, w=4.0)) == ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)
    @test convert(ProbabilityDistribution{Univariate, GaussianMeanPrecision}, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=0.25)) == ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)
    @test convert(ProbabilityDistribution{Multivariate, GaussianMeanPrecision}, ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[8.0], w=mat(4.0))) == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))
    @test convert(ProbabilityDistribution{Multivariate, GaussianMeanPrecision}, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.25))) == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))
    @test convert(ProbabilityDistribution{Multivariate, GaussianMeanPrecision}, ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=1.0, w=2.0)) == ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0))
end


#-------------
# Update rules
#-------------

@testset "SPGaussianMeanPrecisionOutNPP" begin
    @test SPGaussianMeanPrecisionOutNPP <: SumProductRule{GaussianMeanPrecision}
    @test outboundType(SPGaussianMeanPrecisionOutNPP) == Message{GaussianMeanPrecision}
    @test isApplicable(SPGaussianMeanPrecisionOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanPrecisionOutNPP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanPrecisionOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanPrecision, m=1.0, w=2.0)
    @test ruleSPGaussianMeanPrecisionOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0))
end

@testset "SPGaussianMeanPrecisionMPNP" begin
    @test SPGaussianMeanPrecisionMPNP <: SumProductRule{GaussianMeanPrecision}
    @test outboundType(SPGaussianMeanPrecisionMPNP) == Message{GaussianMeanPrecision}
    @test !isApplicable(SPGaussianMeanPrecisionMPNP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanPrecisionMPNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanPrecisionMPNP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanPrecision, m=1.0, w=2.0)
    @test ruleSPGaussianMeanPrecisionMPNP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0))
end

@testset "SPGaussianMeanPrecisionOutNGP" begin
    @test SPGaussianMeanPrecisionOutNGP <: SumProductRule{GaussianMeanPrecision}
    @test outboundType(SPGaussianMeanPrecisionOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanPrecisionOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanPrecisionOutNGP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanPrecisionOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=1.5)
    @test ruleSPGaussianMeanPrecisionOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.5))
end

@testset "SPGaussianMeanPrecisionMGNP" begin
    @test SPGaussianMeanPrecisionMGNP <: SumProductRule{GaussianMeanPrecision}
    @test outboundType(SPGaussianMeanPrecisionMGNP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanPrecisionMGNP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanPrecisionMGNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanPrecisionMGNP(Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=1.5)
    @test ruleSPGaussianMeanPrecisionMGNP(Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.5))
end

@testset "VBGaussianMeanPrecisionM" begin
    @test VBGaussianMeanPrecisionM <: NaiveVariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionM) == Message{GaussianMeanPrecision}
    @test isApplicable(VBGaussianMeanPrecisionM, [ProbabilityDistribution, Nothing, ProbabilityDistribution])
    @test !isApplicable(VBGaussianMeanPrecisionM, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBGaussianMeanPrecisionM(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanPrecision, m=3.0, w=0.5)
    @test ruleVBGaussianMeanPrecisionM(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, ProbabilityDistribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(0.5))
end

@testset "VBGaussianMeanPrecisionW" begin
    @test VBGaussianMeanPrecisionW <: NaiveVariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(VBGaussianMeanPrecisionW, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBGaussianMeanPrecisionW(ProbabilityDistribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing) == Message(Univariate, Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2))
    @test ruleVBGaussianMeanPrecisionW(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing) == Message(MatrixVariate, Wishart, v=mat(1.0/(2.0 + 4.0 + (3.0 - 1.0)^2)), nu=3.0)
end

@testset "VBGaussianMeanPrecisionOut" begin
    @test VBGaussianMeanPrecisionOut <: NaiveVariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionOut) == Message{GaussianMeanPrecision}
    @test isApplicable(VBGaussianMeanPrecisionOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])

    @test ruleVBGaussianMeanPrecisionOut(nothing, ProbabilityDistribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanPrecision, m=3.0, w=0.5)
    @test ruleVBGaussianMeanPrecisionOut(nothing, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), ProbabilityDistribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(0.5))
end

@testset "SVBGaussianMeanPrecisionMGVD" begin
    @test SVBGaussianMeanPrecisionMGVD <: StructuredVariationalRule{GaussianMeanPrecision}
    @test outboundType(SVBGaussianMeanPrecisionMGVD) == Message{GaussianMeanVariance}
    @test isApplicable(SVBGaussianMeanPrecisionMGVD, [Message{Gaussian}, Nothing, ProbabilityDistribution])

    @test ruleSVBGaussianMeanPrecisionMGVD(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanVariance, m=3.0, v=6.0)
    @test ruleSVBGaussianMeanPrecisionMGVD(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, ProbabilityDistribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(6.0))
end

@testset "SVBGaussianMeanPrecisionW" begin
    @test SVBGaussianMeanPrecisionW <: StructuredVariationalRule{GaussianMeanPrecision}
    @test outboundType(SVBGaussianMeanPrecisionW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(SVBGaussianMeanPrecisionW, [ProbabilityDistribution, Nothing])

    @test ruleSVBGaussianMeanPrecisionW(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[2.0, 3.0], v=[5.0 1.0; 1.0 4.0]), nothing) == Message(Univariate, Gamma, a=1.5, b=0.5*(5.0 - 2*1.0 + 4.0 + (3.0 - 2.0)^2))
    @test ruleSVBGaussianMeanPrecisionW(ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.0, 2.0, 3.0, 4.0], v=[5.0 1.0 0.5 0.0; 1.0 4.0 2.0 0.5; 0.5 2.0 3.0 1.0; 0.0 0.5 1.0 2.0]), nothing) == Message(MatrixVariate, Wishart, v=cholinv([5.0 1.0; 1.0 4.0] - [0.5 0.0; 2.0 0.5] - [0.5 2.0; 0.0 0.5] + [3.0 1.0; 1.0 2.0] + ([1.0, 2.0] - [3.0, 4.0])*([1.0, 2.0] - [3.0, 4.0])'), nu=4.0)
end

@testset "SVBGaussianMeanPrecisionOutVGD" begin
    @test SVBGaussianMeanPrecisionOutVGD <: StructuredVariationalRule{GaussianMeanPrecision}
    @test outboundType(SVBGaussianMeanPrecisionOutVGD) == Message{GaussianMeanVariance}
    @test isApplicable(SVBGaussianMeanPrecisionOutVGD, [Nothing, Message{Gaussian}, ProbabilityDistribution])

    @test ruleSVBGaussianMeanPrecisionOutVGD(nothing, Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanVariance, m=3.0, v=6.0)
    @test ruleSVBGaussianMeanPrecisionOutVGD(nothing, Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), ProbabilityDistribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(6.0))
end

@testset "MGaussianMeanPrecisionGGD" begin
    @test MGaussianMeanPrecisionGGD <: MarginalRule{GaussianMeanPrecision}
    @test isApplicable(MGaussianMeanPrecisionGGD, [Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution])

    @test ruleMGaussianMeanPrecisionGGD(Message(Univariate, GaussianMeanPrecision, m=1.0, w=2.0), Message(Univariate, GaussianMeanPrecision, m=3.0, w=4.0), ProbabilityDistribution(Univariate, Gamma, a=1.0, b=2.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
    @test ruleMGaussianMeanPrecisionGGD(Message(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0)), Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0)), ProbabilityDistribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0)) == averageEnergy(GaussianMeanPrecision, ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0), ProbabilityDistribution(Univariate, PointMass, m=0.0), ProbabilityDistribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0)) == differentialEntropy(ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(2.0)))
    @test averageEnergy(GaussianMeanPrecision, ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0), ProbabilityDistribution(Univariate, PointMass, m=0.0), ProbabilityDistribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanPrecision, ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(2.0)), ProbabilityDistribution(Multivariate, PointMass, m=[0.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)))
    @test averageEnergy(GaussianMeanPrecision, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), ProbabilityDistribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanPrecision, ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(2.0)))
end

end #module
