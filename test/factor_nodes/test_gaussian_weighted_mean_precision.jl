module GaussianWeightedMeanPrecisionTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, isProper, unsafeMean, unsafeMode, unsafeVar, unsafeCov, unsafeMeanCov, unsafePrecision, unsafeMeanPrecision, unsafeWeightedMean, unsafeWeightedMeanPrecision
using ForneyLab: SPGaussianWeightedMeanPrecisionOutNPP, VBGaussianWeightedMeanPrecisionOut

@testset "dims" begin
    @test dims(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)) == ()
    @test dims(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=ones(2), w=diageye(2))) == (2,)
end

@testset "vague" begin
    @test vague(GaussianWeightedMeanPrecision) == Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
    @test vague(GaussianWeightedMeanPrecision, (2,)) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(2), w=tiny*eye(2))
end

@testset "isProper" begin
    # Univariate
    @test isProper(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0))
    @test !isProper(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=-1.0))

    # Multivariate
    @test isProper(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0)))
    @test isProper(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=ones(2), w=diageye(2)))
    @test !isProper(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(-1.0)))
end

@testset "==" begin
    # Univariate
    @test Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) == Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0)
    @test Distribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0) == Distribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)

    # Multivariate
    @test Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0)) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0))
    @test Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[0.0], w=mat(1.0)) == Distribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))
end

@testset "unsafe statistics" begin
    # Univariate
    @test unsafeMean(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.5
    @test unsafeMode(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.5
    @test unsafeVar(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.25
    @test unsafeCov(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 0.25
    @test unsafeMeanCov(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == (0.5, 0.25)
    @test unsafePrecision(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 4.0
    @test unsafeMeanPrecision(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == (0.5, 4.0)
    @test unsafeWeightedMean(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == 2.0
    @test unsafeWeightedMeanPrecision(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)) == (2.0, 4.0)

    # Multivariate
    @test unsafeMean(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [0.5]
    @test unsafeMode(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [0.5]
    @test unsafeVar(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [0.25]
    @test unsafeCov(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == mat(0.25)
    @test unsafeMeanCov(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == ([0.5], mat(0.25))
    @test unsafePrecision(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == mat(4.0)
    @test unsafeMeanPrecision(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == ([0.5], mat(4.0))
    @test unsafeWeightedMean(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == [2.0]
    @test unsafeWeightedMeanPrecision(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))) == ([2.0], mat(4.0))
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=2.0), 1.0), -0.5723649429247)
    @test isapprox(logPdf(Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0, 2.0], w=[2.0 0.0; 0.0 2.0]), [1.0, 0.0]), -2.1447298858494)
end

@testset "convert" begin
    @test convert(Distribution{Univariate, GaussianWeightedMeanPrecision}, Distribution(Univariate, GaussianMeanPrecision, m=0.5, w=4.0)) == Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)
    @test convert(Distribution{Univariate, GaussianWeightedMeanPrecision}, Distribution(Univariate, GaussianMeanVariance, m=0.5, v=0.25)) == Distribution(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=4.0)
    @test convert(Distribution{Multivariate, GaussianWeightedMeanPrecision}, Distribution(Multivariate, GaussianMeanPrecision, m=[0.5], w=mat(4.0))) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, GaussianWeightedMeanPrecision}, Distribution(Multivariate, GaussianMeanVariance, m=[0.5], v=mat(0.25))) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(4.0))
    @test convert(Distribution{Multivariate, GaussianWeightedMeanPrecision}, Distribution(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0)) == Distribution(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0))
end

#-------------
# Update rules
#-------------

@testset "SPGaussianWeightedMeanPrecisionOutNPP" begin
    @test SPGaussianWeightedMeanPrecisionOutNPP <: SumProductRule{GaussianWeightedMeanPrecision}
    @test outboundType(SPGaussianWeightedMeanPrecisionOutNPP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPGaussianWeightedMeanPrecisionOutNPP, [Nothing, Message{PointMass}, Message{PointMass}]) 
    @test !isApplicable(SPGaussianWeightedMeanPrecisionOutNPP, [Message{PointMass}, Nothing, Message{PointMass}]) 

    @test ruleSPGaussianWeightedMeanPrecisionOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0)
    @test ruleSPGaussianWeightedMeanPrecisionOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0))
end

@testset "VBGaussianWeightedMeanPrecisionOut" begin
    @test VBGaussianWeightedMeanPrecisionOut <: NaiveVariationalRule{GaussianWeightedMeanPrecision}
    @test outboundType(VBGaussianWeightedMeanPrecisionOut) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(VBGaussianWeightedMeanPrecisionOut, [Nothing, Distribution, Distribution])
    @test !isApplicable(VBGaussianWeightedMeanPrecisionOut, [Distribution, Distribution, Nothing]) 

    @test ruleVBGaussianWeightedMeanPrecisionOut(nothing, Distribution(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Distribution(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=3.0)
    @test ruleVBGaussianWeightedMeanPrecisionOut(nothing, Distribution(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), Distribution(MatrixVariate, PointMass, m=mat(3.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(3.0))
end

end # module