module GaussianCanonicalTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeMeanPrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianCanonicalOutNPP, VBGaussianCanonicalOut

@testset "dims" begin
    @test dims(Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0)) == ()
    @test dims(Distribution(Multivariate, Gaussian{Canonical}, xi=ones(2), w=diageye(2))) == (2,)
end

@testset "vague" begin
    @test vague(Gaussian{Canonical}) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=tiny)
    @test vague(Gaussian{Canonical}, (2,)) == Distribution(Multivariate, Gaussian{Canonical}, xi=zeros(2), w=tiny*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0))
    @test !isProper(Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=-1.0))

    # Multivariate
    @test isProper(Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(1.0)))
    @test isProper(Distribution(Multivariate, Gaussian{Canonical}, xi=ones(2), w=diageye(2)))
    @test !isProper(Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0) == Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0)
    @test Distribution(Univariate, Gaussian{Canonical}, xi=0.0, w=1.0) == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)

    # Multivariate
    @test Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(1.0)) == Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(1.0))
    @test Distribution(Multivariate, Gaussian{Canonical}, xi=[0.0], w=mat(1.0)) == Distribution(Multivariate, Gaussian{Moments}, m=[0.0], v=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == 0.5
    @test unsafeMode(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == 0.5
    @test unsafeVar(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == 0.25
    @test unsafeCov(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == 0.25
    @test unsafeMeanCov(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == (0.5, 0.25)
    @test unsafePrecision(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == 4.0
    @test unsafeMeanPrecision(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == (0.5, 4.0)
    @test unsafeWeightedMean(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == 2.0
    @test unsafeWeightedMeanPrecision(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)) == (2.0, 4.0)

    # Multivariate
    @test unsafeMean(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == [0.5]
    @test unsafeMode(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == [0.5]
    @test unsafeVar(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == [0.25]
    @test unsafeCov(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == mat(0.25)
    @test unsafeMeanCov(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == ([0.5], mat(0.25))
    @test unsafePrecision(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == mat(4.0)
    @test unsafeMeanPrecision(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == ([0.5], mat(4.0))
    @test unsafeWeightedMean(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeWeightedMeanPrecision(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))) == ([2.0], mat(4.0))
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=2.0), 1.0), -0.5723649429247)
    @test isapprox(logPdf(Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0, 2.0], w=[2.0 0.0; 0.0 2.0]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(Distribution{Univariate, Gaussian{Canonical}}, Distribution(Univariate, Gaussian{Precision}, m=0.5, w=4.0)) == Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)
    @test convert(Distribution{Univariate, Gaussian{Canonical}}, Distribution(Univariate, Gaussian{Moments}, m=0.5, v=0.25)) == Distribution(Univariate, Gaussian{Canonical}, xi=2.0, w=4.0)
    @test convert(Distribution{Multivariate, Gaussian{Canonical}}, Distribution(Multivariate, Gaussian{Precision}, m=[0.5], w=mat(4.0))) == Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, Gaussian{Canonical}}, Distribution(Multivariate, Gaussian{Moments}, m=[0.5], v=mat(0.25))) == Distribution(Multivariate, Gaussian{Canonical}, xi=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, Gaussian{Canonical}}, Distribution(Univariate, Gaussian{Canonical}, xi=1.0, w=2.0)) == Distribution(Multivariate, Gaussian{Canonical}, xi=[1.0], w=mat(2.0))
end

#-------------
# Update rules
#-------------

@testset "SPGaussianCanonicalOutNPP" begin
    @test SPGaussianCanonicalOutNPP <: SumProductRule{Gaussian{Canonical}}
    @test outboundType(SPGaussianCanonicalOutNPP) == Message{Gaussian{Canonical}}
    @test isApplicable(SPGaussianCanonicalOutNPP, [Nothing, Message{PointMass}, Message{PointMass}]) 
    @test !isApplicable(SPGaussianCanonicalOutNPP, [Message{PointMass}, Nothing, Message{PointMass}]) 

    @test ruleSPGaussianCanonicalOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian{Canonical}, xi=1.0, w=2.0)
    @test ruleSPGaussianCanonicalOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian{Canonical}, xi=[1.0], w=mat(2.0))
end

@testset "VBGaussianCanonicalOut" begin
    @test VBGaussianCanonicalOut <: NaiveVariationalRule{Gaussian{Canonical}}
    @test outboundType(VBGaussianCanonicalOut) == Message{Gaussian{Canonical}}
    @test isApplicable(VBGaussianCanonicalOut, [Nothing, Distribution, Distribution])
    @test !isApplicable(VBGaussianCanonicalOut, [Distribution, Distribution, Nothing]) 

    @test ruleVBGaussianCanonicalOut(nothing, Distribution(Univariate, Gaussian{Moments}, m=1.0, v=2.0), Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian{Canonical}, xi=1.0, w=3.0)
    @test ruleVBGaussianCanonicalOut(nothing, Distribution(Multivariate, Gaussian{Moments}, m=[1.0], v=mat(2.0)), Distribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, Gaussian{Canonical}, xi=[1.0], w=mat(3.0))
end

end # module