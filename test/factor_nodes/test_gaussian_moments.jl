module GaussianMeanVarianceTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeMeanPrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianMeanVarianceOutNGS, SPGaussianMeanVarianceOutNPP,SPGaussianMeanVarianceMSNP, SPGaussianMeanVarianceMPNP, SPGaussianMeanVarianceOutNGP, SPGaussianMeanVarianceMGNP, SPGaussianMeanVarianceVGGN, SPGaussianMeanVarianceVPGN, SPGaussianMeanVarianceOutNSP, VBGaussianMeanVarianceM, VBGaussianMeanVarianceOut, bootstrap

@testset "dims" begin
    @test dims(Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == ()
    @test dims(Distribution(Multivariate, GaussianMeanVariance, m=ones(2), v=diageye(2))) == (2,)
end

@testset "vague" begin
    @test vague(GaussianMeanVariance) == Distribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
    @test vague(GaussianMeanVariance, (2,)) == Distribution(Multivariate, GaussianMeanVariance, m=zeros(2), v=huge*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0))
    @test !isProper(Distribution(Univariate, GaussianMeanVariance, m=0.0, v=-1.0))

    # Multivariate
    @test isProper(Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)))
    @test isProper(Distribution(Multivariate, GaussianMeanVariance, m=ones(2), v=diageye(2)))
    @test !isProper(Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    @test Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0) == Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)

    # Multivariate
    @test Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)) == Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
    @test Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0)) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 2.0
    @test unsafeMode(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 2.0
    @test unsafeVar(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeCov(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 4.0
    @test unsafeMeanCov(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (2.0, 4.0)
    @test unsafePrecision(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.25
    @test unsafeMeanPrecision(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (2.0, 0.25)
    @test unsafeWeightedMean(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == 0.5
    @test unsafeWeightedMeanPrecision(Distribution(Univariate, GaussianMeanVariance, m=2.0, v=4.0)) == (0.5, 0.25)

    # Multivariate
    @test unsafeMean(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeMode(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeVar(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [4.0]
    @test unsafeCov(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(4.0)
    @test unsafeMeanCov(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([2.0], mat(4.0))
    @test unsafePrecision(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == mat(0.25)
    @test unsafeMeanPrecision(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([2.0], mat(0.25))
    @test unsafeWeightedMean(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == [0.5]
    @test unsafeWeightedMeanPrecision(Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))) == ([0.5], mat(0.25))
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, GaussianMeanVariance, m=1.0, v=0.5), 1.0), -0.5723649429247)
    @test isapprox(logPdf(Distribution(Multivariate, GaussianMeanVariance, m=[1.0, 1.0], v=[0.5 0.0; 0.0 0.5]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(Distribution{Univariate, GaussianMeanVariance}, Distribution(Univariate, GaussianWeightedMeanPrecision, xi=8.0, w=4.0)) == Distribution(Univariate, GaussianMeanVariance, m=2.0, v=0.25)
    @test convert(Distribution{Univariate, GaussianMeanVariance}, Distribution(Univariate, GaussianMeanPrecision, m=2.0, w=4.0)) == Distribution(Univariate, GaussianMeanVariance, m=2.0, v=0.25)
    @test convert(Distribution{Multivariate, GaussianMeanVariance}, Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[8.0], w=mat(4.0))) == Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.25))
    @test convert(Distribution{Multivariate, GaussianMeanVariance}, Distribution(Multivariate, GaussianMeanPrecision, m=[2.0], w=mat(4.0))) == Distribution(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.25))
    @test convert(Distribution{Multivariate, GaussianMeanVariance}, Distribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Distribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end


#-------------
# Update rules
#-------------

@testset "SPGaussianMeanVarianceOutNPP" begin
    @test SPGaussianMeanVarianceOutNPP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutNPP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanVarianceOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutNPP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceMPNP" begin
    @test SPGaussianMeanVarianceMPNP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMPNP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanVarianceMPNP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanVarianceMPNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceMPNP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
    @test ruleSPGaussianMeanVarianceMPNP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMeanVarianceOutNGP" begin
    @test SPGaussianMeanVarianceOutNGP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPGaussianMeanVarianceOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutNGP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMeanVarianceMGNP" begin
    @test SPGaussianMeanVarianceMGNP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMGNP) == Message{GaussianMeanVariance}
    @test !isApplicable(SPGaussianMeanVarianceMGNP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test isApplicable(SPGaussianMeanVarianceMGNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceMGNP(Message(Univariate, GaussianMeanVariance, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleSPGaussianMeanVarianceMGNP(Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMeanVarianceVGGN" begin
    @test SPGaussianMeanVarianceVGGN <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVGGN) == Message{Function}
    @test !isApplicable(SPGaussianMeanVarianceVGGN, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test isApplicable(SPGaussianMeanVarianceVGGN, [Message{Gaussian}, Message{Gaussian}, Nothing])

    msg = ruleSPGaussianMeanVarianceVGGN(Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing)
    @test isa(msg, Message{Function, Univariate})
    @test msg.dist.params[:log_pdf](1.0) == -0.5*log(2.0 + 4.0 + 1.0) - 1/(2*1.0)*(1.0 - 3.0)^2
end

@testset "SPGaussianMeanVarianceVPGN" begin
    @test SPGaussianMeanVarianceVPGN <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceVPGN) == Message{Function}
    @test isApplicable(SPGaussianMeanVarianceVPGN, [Message{PointMass}, Message{Gaussian}, Nothing])

    msg = ruleSPGaussianMeanVarianceVPGN(Message(Univariate, PointMass, m=1.0), Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing)
    @test isa(msg, Message{Function, Univariate})
    @test msg.dist.params[:log_pdf](1.0) == -0.5*log(4.0 + 1.0) - 1/(2*1.0)*(1.0 - 3.0)^2
end

@testset "SPGaussianMeanVarianceOutNSP" begin
    @test SPGaussianMeanVarianceOutNSP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutNSP) == Message{SampleList}
    @test isApplicable(SPGaussianMeanVarianceOutNSP, [Nothing, Message{SampleList}, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceOutNSP, [Message{SampleList}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMeanVarianceOutNSP(nothing, Message(Univariate, SampleList, s=[2.0], w=1.0), Message(Univariate, PointMass, m=0.0)) == Message(Univariate, SampleList, s=[2.0], w=1.0)
    msg = ruleSPGaussianMeanVarianceOutNSP(nothing, Message(Multivariate, SampleList, s=[[2.0]], w=[1.0]), Message(MatrixVariate, PointMass, m=mat(tiny)))
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "SPGaussianMeanVarianceMSNP" begin
    @test SPGaussianMeanVarianceMSNP <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceMSNP) == Message{SampleList}
    @test isApplicable(SPGaussianMeanVarianceMSNP, [Message{SampleList}, Nothing, Message{PointMass}])
    @test !isApplicable(SPGaussianMeanVarianceMSNP, [Message{Gaussian}, Nothing, Message{PointMass}])
end

@testset "SPGaussianMeanVarianceOutNGS" begin
    @test SPGaussianMeanVarianceOutNGS <: SumProductRule{GaussianMeanVariance}
    @test outboundType(SPGaussianMeanVarianceOutNGS) == Message{SampleList}
    @test isApplicable(SPGaussianMeanVarianceOutNGS, [Nothing, Message{Gaussian}, Message{SampleList}])
    @test !isApplicable(SPGaussianMeanVarianceOutNGS, [Message{Gaussian}, Nothing, Message{SampleList}])

    @test ruleSPGaussianMeanVarianceOutNGS(nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=0.0), Message(Univariate, SampleList, s=[0.0], w=[1.0])) == Message(Univariate, SampleList, s=[2.0], w=[1.0])
    msg = ruleSPGaussianMeanVarianceOutNGS(nothing, Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(0.0)), Message(MatrixVariate, SampleList, s=[mat(tiny)], w=[1.0]))
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "VBGaussianMeanVarianceM" begin
    @test VBGaussianMeanVarianceM <: NaiveVariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceM) == Message{GaussianMeanVariance}
    @test isApplicable(VBGaussianMeanVarianceM, [Distribution, Nothing, Distribution])
    @test !isApplicable(VBGaussianMeanVarianceM, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianMeanVarianceM(Distribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing, Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceM(Distribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing, Distribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "VBGaussianMeanVarianceOut" begin
    @test VBGaussianMeanVarianceOut <: NaiveVariationalRule{GaussianMeanVariance}
    @test outboundType(VBGaussianMeanVarianceOut) == Message{GaussianMeanVariance}
    @test isApplicable(VBGaussianMeanVarianceOut, [Nothing, Distribution, Distribution])
    @test !isApplicable(VBGaussianMeanVarianceOut, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianMeanVarianceOut(nothing, Distribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    @test ruleVBGaussianMeanVarianceOut(nothing, Distribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), Distribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0)) == averageEnergy(GaussianMeanVariance, Distribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(Distribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0)) == differentialEntropy(Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(2.0)))
    @test averageEnergy(GaussianMeanVariance, Distribution(Univariate, GaussianMeanVariance, m=0.0, v=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0)) == averageEnergy(GaussianMeanVariance, Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(2.0)), Distribution(Multivariate, PointMass, m=[0.0]), Distribution(MatrixVariate, PointMass, m=mat(2.0)))
    @test averageEnergy(GaussianMeanVariance, Distribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(Univariate, PointMass, m=0.5)) == averageEnergy(GaussianMeanPrecision, Distribution(Multivariate, GaussianMeanVariance, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(Univariate, PointMass, m=2.0))
end

end # module