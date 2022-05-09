module GaussianMomentsTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeMeanPrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianMomentsOutNGS, SPGaussianMomentsOutNPP,SPGaussianMomentsMSNP, SPGaussianMomentsMPNP, SPGaussianMomentsOutNGP, SPGaussianMomentsMGNP, SPGaussianMomentsVGGN, SPGaussianMomentsVPGN, SPGaussianMomentsOutNSP, VBGaussianMomentsM, VBGaussianMomentsOut, bootstrap

@testset "default Gaussian{Moments} node definition" begin
    fg = FactorGraph()
    x = Variable()
    nd = Gaussian(x, constant(0.0), constant(1.0))
    @test isa(nd, Gaussian{Moments})
end

@testset "default Gaussian{Moments} Distribution definitions" begin
    @test Distribution(Univariate, Gaussian; m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test Distribution(Gaussian; m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test Distribution(Multivariate, Gaussian; m=[0.0], v=mat(1.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0))
end

@testset "dims" begin
    @test dims(Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)) == ()
    @test dims(Distribution(Multivariate, Gaussian{Moments}, m=ones(2), v=diageye(2))) == (2,)
end

@testset "vague" begin
    @test vague(Gaussian{Moments}) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=huge)
    @test vague(Gaussian{Moments}, (2,)) == Distribution(Multivariate, Gaussian{Moments}, m=zeros(2), v=huge*eye(2))
    @test vague(Gaussian) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=huge)
    @test vague(Gaussian, (2,)) == Distribution(Multivariate, Gaussian{Moments}, m=zeros(2), v=huge*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0))
    @test !isProper(Distribution(Univariate, Gaussian{Moments}, m=0.0, v=-1.0))

    # Multivariate
    @test isProper(Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0)))
    @test isProper(Distribution(Multivariate, Gaussian{Moments}, m=ones(2), v=diageye(2)))
    @test !isProper(Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0)

    # Multivariate
    @test Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0))
    @test Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0)) == Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == 2.0
    @test unsafeMode(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == 2.0
    @test unsafeVar(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == 4.0
    @test unsafeCov(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == 4.0
    @test unsafeMeanCov(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == (2.0, 4.0)
    @test unsafePrecision(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == 0.25
    @test unsafeMeanPrecision(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == (2.0, 0.25)
    @test unsafeWeightedMean(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == 0.5
    @test unsafeWeightedMeanPrecision(Distribution(Univariate, Gaussian{Moments}, m=2.0, v=4.0)) == (0.5, 0.25)

    # Multivariate
    @test unsafeMean(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeMode(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == [2.0]
    @test unsafeVar(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == [4.0]
    @test unsafeCov(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == mat(4.0)
    @test unsafeMeanCov(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == ([2.0], mat(4.0))
    @test unsafePrecision(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == mat(0.25)
    @test unsafeMeanPrecision(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == ([2.0], mat(0.25))
    @test unsafeWeightedMean(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == [0.5]
    @test unsafeWeightedMeanPrecision(Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(4.0))) == ([0.5], mat(0.25))
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, Gaussian{Moments}, m=1.0, v=0.5), 1.0), -0.5723649429247)
    @test isapprox(logPdf(Distribution(Multivariate, Gaussian{Moments}, m=[1.0, 1.0], v=[0.5 0.0; 0.0 0.5]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(Distribution{Univariate, Gaussian{Moments}}, Distribution(Univariate, Gaussian{Canonical}, xi=8.0, w=4.0)) == Distribution(Univariate, Gaussian{Moments}, m=2.0, v=0.25)
    @test convert(Distribution{Univariate, Gaussian{Moments}}, Distribution(Univariate, Gaussian{Precision}, m=2.0, w=4.0)) == Distribution(Univariate, Gaussian{Moments}, m=2.0, v=0.25)
    @test convert(Distribution{Multivariate, Gaussian{Moments}}, Distribution(Multivariate, Gaussian{Canonical}, xi=[8.0], w=mat(4.0))) == Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(0.25))
    @test convert(Distribution{Multivariate, Gaussian{Moments}}, Distribution(Multivariate, Gaussian{Precision}, m=[2.0], w=mat(4.0))) == Distribution(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(0.25))
    @test convert(Distribution{Multivariate, Gaussian{Moments}}, Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0))
end


#-------------
# Update rules
#-------------

@testset "SPGaussianMomentsOutNPP" begin
    @test SPGaussianMomentsOutNPP <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsOutNPP) == Message{Gaussian{Moments}}
    @test isApplicable(SPGaussianMomentsOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPGaussianMomentsOutNPP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMomentsOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=2.0)
    @test ruleSPGaussianMomentsOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMomentsMPNP" begin
    @test SPGaussianMomentsMPNP <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsMPNP) == Message{Gaussian{Moments}}
    @test !isApplicable(SPGaussianMomentsMPNP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test isApplicable(SPGaussianMomentsMPNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMomentsMPNP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=2.0)
    @test ruleSPGaussianMomentsMPNP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0))
end

@testset "SPGaussianMomentsOutNGP" begin
    @test SPGaussianMomentsOutNGP <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsOutNGP) == Message{Gaussian{Moments}}
    @test isApplicable(SPGaussianMomentsOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test !isApplicable(SPGaussianMomentsOutNGP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMomentsOutNGP(nothing, Message(Univariate, Gaussian{Moments}, m=1.0, v=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=3.0)
    @test ruleSPGaussianMomentsOutNGP(nothing, Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(1.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMomentsMGNP" begin
    @test SPGaussianMomentsMGNP <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsMGNP) == Message{Gaussian{Moments}}
    @test !isApplicable(SPGaussianMomentsMGNP, [Nothing, Message{Gaussian}, Message{PointMass}])
    @test isApplicable(SPGaussianMomentsMGNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMomentsMGNP(Message(Univariate, Gaussian{Moments}, m=1.0, v=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=3.0)
    @test ruleSPGaussianMomentsMGNP(Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(1.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(3.0))
end

@testset "SPGaussianMomentsVGGN" begin
    @test SPGaussianMomentsVGGN <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsVGGN) == Message{Function}
    @test !isApplicable(SPGaussianMomentsVGGN, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test isApplicable(SPGaussianMomentsVGGN, [Message{Gaussian}, Message{Gaussian}, Nothing])

    msg = ruleSPGaussianMomentsVGGN(Message(Univariate, Gaussian{Moments}, m=1.0, v=2.0), Message(Univariate, Gaussian{Moments}, m=3.0, v=4.0), nothing)
    @test isa(msg, Message{Function, Univariate})
    @test msg.dist.params[:log_pdf](1.0) == -0.5*log(2.0 + 4.0 + 1.0) - 1/(2*(2.0 + 4.0 + 1.0))*(1.0 - 3.0)^2
end

@testset "SPGaussianMomentsVPGN" begin
    @test SPGaussianMomentsVPGN <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsVPGN) == Message{Function}
    @test isApplicable(SPGaussianMomentsVPGN, [Message{PointMass}, Message{Gaussian}, Nothing])

    msg = ruleSPGaussianMomentsVPGN(Message(Univariate, PointMass, m=1.0), Message(Univariate, Gaussian{Moments}, m=3.0, v=4.0), nothing)
    @test isa(msg, Message{Function, Univariate})
    @test msg.dist.params[:log_pdf](1.0) == -0.5*log(4.0 + 1.0) - 1/(2*(4.0 + 1.0))*(1.0 - 3.0)^2
end

@testset "SPGaussianMomentsOutNSP" begin
    @test SPGaussianMomentsOutNSP <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsOutNSP) == Message{SampleList}
    @test isApplicable(SPGaussianMomentsOutNSP, [Nothing, Message{SampleList}, Message{PointMass}])
    @test !isApplicable(SPGaussianMomentsOutNSP, [Message{SampleList}, Nothing, Message{PointMass}])

    @test ruleSPGaussianMomentsOutNSP(nothing, Message(Univariate, SampleList, s=[2.0], w=1.0), Message(Univariate, PointMass, m=0.0)) == Message(Univariate, SampleList, s=[2.0], w=1.0)
    msg = ruleSPGaussianMomentsOutNSP(nothing, Message(Multivariate, SampleList, s=[[2.0]], w=[1.0]), Message(MatrixVariate, PointMass, m=mat(tiny)))
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "SPGaussianMomentsMSNP" begin
    @test SPGaussianMomentsMSNP <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsMSNP) == Message{SampleList}
    @test isApplicable(SPGaussianMomentsMSNP, [Message{SampleList}, Nothing, Message{PointMass}])
    @test !isApplicable(SPGaussianMomentsMSNP, [Message{Gaussian}, Nothing, Message{PointMass}])
end

@testset "SPGaussianMomentsOutNGS" begin
    @test SPGaussianMomentsOutNGS <: SumProductRule{Gaussian{Moments}}
    @test outboundType(SPGaussianMomentsOutNGS) == Message{SampleList}
    @test isApplicable(SPGaussianMomentsOutNGS, [Nothing, Message{Gaussian}, Message{SampleList}])
    @test !isApplicable(SPGaussianMomentsOutNGS, [Message{Gaussian}, Nothing, Message{SampleList}])

    @test ruleSPGaussianMomentsOutNGS(nothing, Message(Univariate, Gaussian{Moments}, m=2.0, v=0.0), Message(Univariate, SampleList, s=[0.0], w=[1.0])) == Message(Univariate, SampleList, s=[2.0], w=[1.0])
    msg = ruleSPGaussianMomentsOutNGS(nothing, Message(Multivariate, Gaussian{Moments}, m=[2.0], v=mat(0.0)), Message(MatrixVariate, SampleList, s=[mat(tiny)], w=[1.0]))
    @test isapprox(msg.dist.params[:s][1][1], 2.0, atol=1e-4)
    @test msg.dist.params[:w] == [1.0]
end

@testset "VBGaussianMomentsM" begin
    @test VBGaussianMomentsM <: NaiveVariationalRule{Gaussian{Moments}}
    @test outboundType(VBGaussianMomentsM) == Message{Gaussian{Moments}}
    @test isApplicable(VBGaussianMomentsM, [Distribution, Nothing, Distribution])
    @test !isApplicable(VBGaussianMomentsM, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianMomentsM(Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0), nothing, Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=3.0)
    @test ruleVBGaussianMomentsM(Distribution(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0)), nothing, Distribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(3.0))
end

@testset "VBGaussianMomentsOut" begin
    @test VBGaussianMomentsOut <: NaiveVariationalRule{Gaussian{Moments}}
    @test outboundType(VBGaussianMomentsOut) == Message{Gaussian{Moments}}
    @test isApplicable(VBGaussianMomentsOut, [Nothing, Distribution, Distribution])
    @test !isApplicable(VBGaussianMomentsOut, [Distribution, Distribution, Nothing])

    @test ruleVBGaussianMomentsOut(nothing, Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0), Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian{Moments}, m=1.0, v=3.0)
    @test ruleVBGaussianMomentsOut(nothing, Distribution(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0)), Distribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, Gaussian{Moments}, m=0.0, v=2.0)) == averageEnergy(Gaussian{Moments}, Distribution(Univariate, Gaussian{Moments}, m=0.0, v=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(Distribution(Univariate, Gaussian{Moments}, m=0.0, v=2.0)) == differentialEntropy(Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(2.0)))
    @test averageEnergy(Gaussian{Moments}, Distribution(Univariate, Gaussian{Moments}, m=0.0, v=2.0), Distribution(Univariate, PointMass, m=0.0), Distribution(Univariate, PointMass, m=2.0)) == averageEnergy(Gaussian{Moments}, Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(2.0)), Distribution(Multivariate, PointMass, m=[0.0]), Distribution(MatrixVariate, PointMass, m=mat(2.0)))
    @test averageEnergy(Gaussian{Moments}, Distribution(Multivariate, Gaussian{Moments}, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(Univariate, PointMass, m=0.5)) == averageEnergy(Gaussian{Precision}, Distribution(Multivariate, Gaussian{Moments}, m=[0.0, 1.0], v=[3.0 1.0; 1.0 2.0]), Distribution(Univariate, PointMass, m=2.0))
end

end # module