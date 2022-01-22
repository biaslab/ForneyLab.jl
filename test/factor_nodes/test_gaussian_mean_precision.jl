module GaussianMeanPrecisionTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeMeanPrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianMeanPrecisionOutNPP, SPGaussianMeanPrecisionMPNP, SPGaussianMeanPrecisionOutNGP, SPGaussianMeanPrecisionMGNP, VBGaussianMeanPrecisionOut, VBGaussianMeanPrecisionM, VBGaussianMeanPrecisionW, SVBGaussianMeanPrecisionOutVGD, SVBGaussianMeanPrecisionMGVD, SVBGaussianMeanPrecisionW, MGaussianMeanPrecisionGGD, MGaussianMeanPrecisionGGN

@testset "dims" begin
    @test dims(Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0)) == ()
    @test dims(Distribution(Multivariate, GaussianMeanPrecision, m=ones(2), w=diageye(2))) == (2,)
end

@testset "vague" begin
    @test vague(GaussianMeanPrecision) == Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=tiny)
    @test vague(GaussianMeanPrecision, (2,)) == Distribution(Multivariate, GaussianMeanPrecision, m=zeros(2), w=tiny*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0))
    @test !isProper(Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=-1.0))

    # Multivariate
    @test isProper(Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0)))
    @test isProper(Distribution(Multivariate, GaussianMeanPrecision, m=ones(2), w=diageye(2)))
    @test !isProper(Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0) == Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0)
    @test Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0) == Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)

    # Multivariate
    @test Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0)) == Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0))
    @test Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0)) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 2.0
    @test unsafeMode(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 2.0
    @test unsafeVar(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 0.25
    @test unsafeCov(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 0.25
    @test unsafeMeanCov(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == (2.0, 0.25)
    @test unsafePrecision(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 4.0
    @test unsafeMeanPrecision(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == (2.0, 4.0)
    @test unsafeWeightedMean(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == 8.0
    @test unsafeWeightedMeanPrecision(Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == (8.0, 4.0)

    # Multivariate
    @test unsafeMean(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeMode(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeVar(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [0.25]
    @test unsafeCov(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == mat(0.25)
    @test unsafeMeanCov(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == ([2.0], mat(0.25))
    @test unsafePrecision(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == mat(4.0)
    @test unsafeMeanPrecision(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == ([2.0], mat(4.0))
    @test unsafeWeightedMean(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == [8.0]
    @test unsafeWeightedMeanPrecision(Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == ([8.0], mat(4.0))
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, GaussianMeanPrecision, m=1.0, w=2.0), 1.0), -0.5723649429247)
    @test isapprox(logPdf(Distribution(Multivariate, GaussianMeanPrecision, m=[1.0, 1.0], w=[2.0 0.0; 0.0 2.0]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(Distribution{Univariate, GaussianMeanPrecision}, Distribution(Univariate, GaussianWeightedMeanPrecision, xi=8.0, w=4.0)) == Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)
    @test convert(Distribution{Univariate, GaussianMeanPrecision}, Distribution(Univariate, GaussianMeanVariance, m=2.0, v=0.25)) == Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)
    @test convert(Distribution{Multivariate, GaussianMeanPrecision}, Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[8.0], w=mat(4.0))) == Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, GaussianMeanPrecision}, Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.25))) == Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, GaussianMeanPrecision}, Distribution(Univariate, GaussianMeanPrecision, m=1.0, w=2.0)) == Distribution(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0))
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
    @test isApplicable(VBGaussianMeanPrecisionM, [Distribution, Nothing, Distribution])
    @test !isApplicable(VBGaussianMeanPrecisionM, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianMeanPrecisionM(Distribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanPrecision, m=3.0, w=0.5)
    @test ruleVBGaussianMeanPrecisionM(Distribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(0.5))
end

@testset "VBGaussianMeanPrecisionW" begin
    @test VBGaussianMeanPrecisionW <: NaiveVariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(VBGaussianMeanPrecisionW, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianMeanPrecisionW(Distribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Distribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing) == Message(Univariate, Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2))
    @test ruleVBGaussianMeanPrecisionW(Distribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Distribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing) == Message(MatrixVariate, Wishart, v=mat(1.0/(2.0 + 4.0 + (3.0 - 1.0)^2)), nu=3.0)
end

@testset "VBGaussianMeanPrecisionOut" begin
    @test VBGaussianMeanPrecisionOut <: NaiveVariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecisionOut) == Message{GaussianMeanPrecision}
    @test isApplicable(VBGaussianMeanPrecisionOut, [Nothing, Distribution, Distribution])

    @test ruleVBGaussianMeanPrecisionOut(nothing, Distribution(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanPrecision, m=3.0, w=0.5)
    @test ruleVBGaussianMeanPrecisionOut(nothing, Distribution(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(0.5))
end

@testset "SVBGaussianMeanPrecisionMGVD" begin
    @test SVBGaussianMeanPrecisionMGVD <: StructuredVariationalRule{GaussianMeanPrecision}
    @test outboundType(SVBGaussianMeanPrecisionMGVD) == Message{GaussianMeanVariance}
    @test isApplicable(SVBGaussianMeanPrecisionMGVD, [Message{Gaussian}, Nothing, Distribution])

    @test ruleSVBGaussianMeanPrecisionMGVD(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanVariance, m=3.0, v=6.0)
    @test ruleSVBGaussianMeanPrecisionMGVD(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(6.0))
end

@testset "SVBGaussianMeanPrecisionW" begin
    @test SVBGaussianMeanPrecisionW <: StructuredVariationalRule{GaussianMeanPrecision}
    @test outboundType(SVBGaussianMeanPrecisionW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(SVBGaussianMeanPrecisionW, [Distribution, Nothing])

    @test ruleSVBGaussianMeanPrecisionW(Distribution(Multivariate, GaussianMeanVariance, m=[2.0, 3.0], v=[5.0 1.0; 1.0 4.0]), nothing) == Message(Univariate, Gamma, a=1.5, b=0.5*(5.0 - 2*1.0 + 4.0 + (3.0 - 2.0)^2))
    @test ruleSVBGaussianMeanPrecisionW(Distribution(Multivariate, GaussianMeanVariance, m=[1.0, 2.0, 3.0, 4.0], v=[5.0 1.0 0.5 0.0; 1.0 4.0 2.0 0.5; 0.5 2.0 3.0 1.0; 0.0 0.5 1.0 2.0]), nothing) == Message(MatrixVariate, Wishart, v=cholinv([5.0 1.0; 1.0 4.0] - [0.5 0.0; 2.0 0.5] - [0.5 2.0; 0.0 0.5] + [3.0 1.0; 1.0 2.0] + ([1.0, 2.0] - [3.0, 4.0])*([1.0, 2.0] - [3.0, 4.0])'), nu=4.0)
end

@testset "SVBGaussianMeanPrecisionOutVGD" begin
    @test SVBGaussianMeanPrecisionOutVGD <: StructuredVariationalRule{GaussianMeanPrecision}
    @test outboundType(SVBGaussianMeanPrecisionOutVGD) == Message{GaussianMeanVariance}
    @test isApplicable(SVBGaussianMeanPrecisionOutVGD, [Nothing, Message{Gaussian}, Distribution])

    @test ruleSVBGaussianMeanPrecisionOutVGD(nothing, Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, GaussianMeanVariance, m=3.0, v=6.0)
    @test ruleSVBGaussianMeanPrecisionOutVGD(nothing, Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(6.0))
end

@testset "MGaussianMeanPrecisionGGD" begin
    @test MGaussianMeanPrecisionGGD <: MarginalRule{GaussianMeanPrecision}
    @test isApplicable(MGaussianMeanPrecisionGGD, [Message{Gaussian}, Message{Gaussian}, Distribution])
    @test !isApplicable(MGaussianMeanPrecisionGGD, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleMGaussianMeanPrecisionGGD(Message(Univariate, GaussianMeanPrecision, m=1.0, w=2.0), Message(Univariate, GaussianMeanPrecision, m=3.0, w=4.0), Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Distribution(Multivariate, GaussianMeanVariance, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
    @test ruleMGaussianMeanPrecisionGGD(Message(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0)), Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0)), Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Distribution(Multivariate, GaussianMeanVariance, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
end

@testset "MGaussianMeanPrecisionGGN" begin
    @test MGaussianMeanPrecisionGGN <: MarginalRule{GaussianMeanPrecision}
    @test isApplicable(MGaussianMeanPrecisionGGN, [Message{Gaussian}, Message{Gaussian}, Nothing])
    @test !isApplicable(MGaussianMeanPrecisionGGN, [Message{Gaussian}, Message{Gaussian}, Distribution])

    @test ruleMGaussianMeanPrecisionGGN(Message(Univariate, GaussianMeanPrecision, m=1.0, w=2.0), Message(Univariate, GaussianMeanPrecision, m=3.0, w=4.0), Message(Univariate, PointMass, m=0.5)) == Distribution(Multivariate, GaussianMeanVariance, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
    @test ruleMGaussianMeanPrecisionGGN(Message(Multivariate, GaussianMeanPrecision, m=[1.0], w=mat(2.0)), Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0)), Message(MatrixVariate, PointMass, m=mat(0.5))) == Distribution(Multivariate, GaussianMeanVariance, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0)) == averageEnergy(GaussianMeanPrecision, Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0)) == differentialEntropy(Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(2.0)))
    @test averageEnergy(GaussianMeanPrecision, Distribution(Univariate, GaussianMeanPrecision, m=0.0, w=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanPrecision, Distribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(2.0)), Distribution(Multivariate, PointMass, m=[0.0]), Distribution(MatrixVariate, PointMass, m=mat(2.0)))
    @test averageEnergy(GaussianMeanPrecision, Distribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanPrecision, Distribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(MatrixVariate, PointMass, m=mat(2.0)))
end

end #module
