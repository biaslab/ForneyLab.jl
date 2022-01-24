module Gaussian{Precision}Test

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeMeanPrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianPrecisionOutNPP, SPGaussianPrecisionMPNP, SPGaussianPrecisionOutNGP, SPGaussianPrecisionMGNP, VBGaussianPrecisionOut, VBGaussianPrecisionM, VBGaussianPrecisionW, SVBGaussianPrecisionOutVGD, SVBGaussianPrecisionMGVD, SVBGaussianPrecisionW, MGaussian{Precision}GGD, MGaussian{Precision}GGN

@testset "dims" begin
    @test dims(Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0)) == ()
    @test dims(Distribution(Multivariate, Gaussian{Precision}, m=ones(2), w=diageye(2))) == (2,)
end

@testset "vague" begin
    @test vague(Gaussian{Precision}) == Distribution(Univariate, Gaussian{Precision}, m=0.0, w=tiny)
    @test vague(Gaussian{Precision}, (2,)) == Distribution(Multivariate, Gaussian{Precision}, m=zeros(2), w=tiny*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0))
    @test !isProper(Distribution(Univariate, Gaussian{Precision}, m=0.0, w=-1.0))

    # Multivariate
    @test isProper(Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(1.0)))
    @test isProper(Distribution(Multivariate, Gaussian{Precision}, m=ones(2), w=diageye(2)))
    @test !isProper(Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0) == Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0)
    @test Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0)

    # Multivariate
    @test Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(1.0)) == Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(1.0))
    @test Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(1.0)) == Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == 2.0
    @test unsafeMode(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == 2.0
    @test unsafeVar(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == 0.25
    @test unsafeCov(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == 0.25
    @test unsafeMeanCov(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == (2.0, 0.25)
    @test unsafePrecision(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == 4.0
    @test unsafeMeanPrecision(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == (2.0, 4.0)
    @test unsafeWeightedMean(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == 8.0
    @test unsafeWeightedMeanPrecision(Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == (8.0, 4.0)

    # Multivariate
    @test unsafeMean(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeMode(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeVar(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == [0.25]
    @test unsafeCov(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == mat(0.25)
    @test unsafeMeanCov(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == ([2.0], mat(0.25))
    @test unsafePrecision(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == mat(4.0)
    @test unsafeMeanPrecision(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == ([2.0], mat(4.0))
    @test unsafeWeightedMean(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == [8.0]
    @test unsafeWeightedMeanPrecision(Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == ([8.0], mat(4.0))
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, Gaussian{Precision}, m=1.0, w=2.0), 1.0), -0.5723649429247)
    @test isapprox(logPdf(Distribution(Multivariate, Gaussian{Precision}, m=[1.0, 1.0], w=[2.0 0.0; 0.0 2.0]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(Distribution{Univariate, Gaussian{Precision}}, Distribution(Univariate, Gaussian{Canonical}, xi=8.0, w=4.0)) == Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)
    @test convert(Distribution{Univariate, Gaussian{Precision}}, Distribution(Univariate, Gaussian{Moments}, m=2.0, v=0.25)) == Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)
    @test convert(Distribution{Multivariate, Gaussian{Precision}}, Distribution(Multivariate, Gaussian{Canonical}, xi=[8.0], w=mat(4.0))) == Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, Gaussian{Precision}}, Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(0.25))) == Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, Gaussian{Precision}}, Distribution(Univariate, Gaussian{Precision}, m=1.0, w=2.0)) == Distribution(Multivariate, Gaussian{Precision}, m=[1.0], w=mat(2.0))
end


#-------------
# Update rules
#-------------

@testset "SPGaussianPrecisionOutNPP" begin
    @test SPGaussianPrecisionOutNPP <: SumProductRule{Gaussian{Precision}}
    @test outboundType(SPGaussianPrecisionOutNPP) == Message{Gaussian{Precision}}
    @test isApplicable(SPGaussianPrecisionOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPGaussianPrecisionOutNPP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianPrecisionOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Precision}, m=1.0, w=2.0)
    @test ruleSPGaussianPrecisionOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Precision}, m=[1.0], w=mat(2.0))
end

@testset "SPGaussianPrecisionMPNP" begin
    @test SPGaussianPrecisionMPNP <: SumProductRule{Gaussian{Precision}}
    @test outboundType(SPGaussianPrecisionMPNP) == Message{Gaussian{Precision}}
    @test !isApplicable(SPGaussianPrecisionMPNP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPGaussianPrecisionMPNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianPrecisionMPNP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Precision}, m=1.0, w=2.0)
    @test ruleSPGaussianPrecisionMPNP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Precision}, m=[1.0], w=mat(2.0))
end

@testset "SPGaussianPrecisionOutNGP" begin
    @test SPGaussianPrecisionOutNGP <: SumProductRule{Gaussian{Precision}}
    @test outboundType(SPGaussianPrecisionOutNGP) == Message{Gaussian{Moments}}
    @test isApplicable(SPGaussianPrecisionOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test !isApplicable(SPGaussianPrecisionOutNGP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianPrecisionOutNGP(nothing, Message(Univariate, Gaussian{Moments}, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=1.5)
    @test ruleSPGaussianPrecisionOutNGP(nothing, Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(1.5))
end

@testset "SPGaussianPrecisionMGNP" begin
    @test SPGaussianPrecisionMGNP <: SumProductRule{Gaussian{Precision}}
    @test outboundType(SPGaussianPrecisionMGNP) == Message{Gaussian{Moments}}
    @test !isApplicable(SPGaussianPrecisionMGNP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test isApplicable(SPGaussianPrecisionMGNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianPrecisionMGNP(Message(Univariate, Gaussian{Moments}, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=1.5)
    @test ruleSPGaussianPrecisionMGNP(Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(1.5))
end

@testset "VBGaussianPrecisionM" begin
    @test VBGaussianPrecisionM <: NaiveVariationalRule{Gaussian{Precision}}
    @test outboundType(VBGaussianPrecisionM) == Message{Gaussian{Precision}}
    @test isApplicable(VBGaussianPrecisionM, [Distribution, Nothing, Distribution])
    @test !isApplicable(VBGaussianPrecisionM, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianPrecisionM(Distribution(Univariate, Gaussian{Moments}, m=3.0, v=4.0), nothing, Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, Gaussian{Precision}, m=3.0, w=0.5)
    @test ruleVBGaussianPrecisionM(Distribution(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(4.0)), nothing, Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, Gaussian{Precision}, m=[3.0], w=mat(0.5))
end

@testset "VBGaussianPrecisionW" begin
    @test VBGaussianPrecisionW <: NaiveVariationalRule{Gaussian{Precision}}
    @test outboundType(VBGaussianPrecisionW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(VBGaussianPrecisionW, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianPrecisionW(Distribution(Univariate, Gaussian{Moments}, m=3.0, v=4.0), Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0), nothing) == Message(Univariate, Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2))
    @test ruleVBGaussianPrecisionW(Distribution(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(4.0)), Distribution(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0)), nothing) == Message(MatrixVariate, Wishart, v=mat(1.0/(2.0 + 4.0 + (3.0 - 1.0)^2)), nu=3.0)
end

@testset "VBGaussianPrecisionOut" begin
    @test VBGaussianPrecisionOut <: NaiveVariationalRule{Gaussian{Precision}}
    @test outboundType(VBGaussianPrecisionOut) == Message{Gaussian{Precision}}
    @test isApplicable(VBGaussianPrecisionOut, [Nothing, Distribution, Distribution])

    @test ruleVBGaussianPrecisionOut(nothing, Distribution(Univariate, Gaussian{Moments}, m=3.0, v=4.0), Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, Gaussian{Precision}, m=3.0, w=0.5)
    @test ruleVBGaussianPrecisionOut(nothing, Distribution(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(4.0)), Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, Gaussian{Precision}, m=[3.0], w=mat(0.5))
end

@testset "SVBGaussianPrecisionMGVD" begin
    @test SVBGaussianPrecisionMGVD <: StructuredVariationalRule{Gaussian{Precision}}
    @test outboundType(SVBGaussianPrecisionMGVD) == Message{Gaussian{Moments}}
    @test isApplicable(SVBGaussianPrecisionMGVD, [Message{Gaussian}, Nothing, Distribution])

    @test ruleSVBGaussianPrecisionMGVD(Message(Univariate, Gaussian{Moments}, m=3.0, v=4.0), nothing, Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, Gaussian{Moments}, m=3.0, v=6.0)
    @test ruleSVBGaussianPrecisionMGVD(Message(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(4.0)), nothing, Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(6.0))
end

@testset "SVBGaussianPrecisionW" begin
    @test SVBGaussianPrecisionW <: StructuredVariationalRule{Gaussian{Precision}}
    @test outboundType(SVBGaussianPrecisionW) == Message{Union{Gamma, Wishart}}
    @test isApplicable(SVBGaussianPrecisionW, [Distribution, Nothing])

    @test ruleSVBGaussianPrecisionW(Distribution(Multivariate, Gaussian{Moments}, m=[2.0, 3.0], v=[5.0 1.0; 1.0 4.0]), nothing) == Message(Univariate, Gamma, a=1.5, b=0.5*(5.0 - 2*1.0 + 4.0 + (3.0 - 2.0)^2))
    @test ruleSVBGaussianPrecisionW(Distribution(Multivariate, Gaussian{Moments}, m=[1.0, 2.0, 3.0, 4.0], v=[5.0 1.0 0.5 0.0; 1.0 4.0 2.0 0.5; 0.5 2.0 3.0 1.0; 0.0 0.5 1.0 2.0]), nothing) == Message(MatrixVariate, Wishart, v=cholinv([5.0 1.0; 1.0 4.0] - [0.5 0.0; 2.0 0.5] - [0.5 2.0; 0.0 0.5] + [3.0 1.0; 1.0 2.0] + ([1.0, 2.0] - [3.0, 4.0])*([1.0, 2.0] - [3.0, 4.0])'), nu=4.0)
end

@testset "SVBGaussianPrecisionOutVGD" begin
    @test SVBGaussianPrecisionOutVGD <: StructuredVariationalRule{Gaussian{Precision}}
    @test outboundType(SVBGaussianPrecisionOutVGD) == Message{Gaussian{Moments}}
    @test isApplicable(SVBGaussianPrecisionOutVGD, [Nothing, Message{Gaussian}, Distribution])

    @test ruleSVBGaussianPrecisionOutVGD(nothing, Message(Univariate, Gaussian{Moments}, m=3.0, v=4.0), Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Message(Univariate, Gaussian{Moments}, m=3.0, v=6.0)
    @test ruleSVBGaussianPrecisionOutVGD(nothing, Message(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(4.0)), Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Message(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(6.0))
end

@testset "MGaussian{Precision}GGD" begin
    @test MGaussian{Precision}GGD <: MarginalRule{Gaussian{Precision}}
    @test isApplicable(MGaussian{Precision}GGD, [Message{Gaussian}, Message{Gaussian}, Distribution])
    @test !isApplicable(MGaussian{Precision}GGD, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleMGaussian{Precision}GGD(Message(Univariate, Gaussian{Precision}, m=1.0, w=2.0), Message(Univariate, Gaussian{Precision}, m=3.0, w=4.0), Distribution(Univariate, Gamma, a=1.0, b=2.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
    @test ruleMGaussian{Precision}GGD(Message(Multivariate, Gaussian{Precision}, m=[1.0], w=mat(2.0)), Message(Multivariate, Gaussian{Precision}, m=[3.0], w=mat(4.0)), Distribution(MatrixVariate, Wishart, v=mat(0.25), nu=2.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
end

@testset "MGaussian{Precision}GGN" begin
    @test MGaussian{Precision}GGN <: MarginalRule{Gaussian{Precision}}
    @test isApplicable(MGaussian{Precision}GGN, [Message{Gaussian}, Message{Gaussian}, Nothing])
    @test !isApplicable(MGaussian{Precision}GGN, [Message{Gaussian}, Message{Gaussian}, Distribution])

    @test ruleMGaussian{Precision}GGN(Message(Univariate, Gaussian{Precision}, m=1.0, w=2.0), Message(Univariate, Gaussian{Precision}, m=3.0, w=4.0), Message(Univariate, PointMass, m=0.5)) == Distribution(Multivariate, Gaussian{Moments}, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
    @test ruleMGaussian{Precision}GGN(Message(Multivariate, Gaussian{Precision}, m=[1.0], w=mat(2.0)), Message(Multivariate, Gaussian{Precision}, m=[3.0], w=mat(4.0)), Message(MatrixVariate, PointMass, m=mat(0.5))) == Distribution(Multivariate, Gaussian{Moments}, m=[1.3636363636363638, 2.8181818181818175], v=[0.4090909090909091 0.04545454545454545; 0.04545454545454545 0.22727272727272724])
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, Gaussian{Precision}, m=0.0, w=2.0)) == averageEnergy(Gaussian{Precision}, Distribution(Univariate, Gaussian{Precision}, m=0.0, w=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(Distribution(Univariate, Gaussian{Precision}, m=0.0, w=2.0)) == differentialEntropy(Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(2.0)))
    @test averageEnergy(Gaussian{Precision}, Distribution(Univariate, Gaussian{Precision}, m=0.0, w=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0)) == averageEnergy(Gaussian{Precision}, Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(2.0)), Distribution(Multivariate, PointMass, m=[0.0]), Distribution(MatrixVariate, PointMass, m=mat(2.0)))
    @test averageEnergy(Gaussian{Precision}, Distribution(Multivariate, Gaussian{Moments}, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(Univariate, PointMass, m=2.0)) == averageEnergy(Gaussian{Precision}, Distribution(Multivariate, Gaussian{Moments}, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(MatrixVariate, PointMass, m=mat(2.0)))
end

end #module
