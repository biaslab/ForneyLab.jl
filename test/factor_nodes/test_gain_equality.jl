module GainEqualityTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPGainEqualityOutVGGP, SPGainEqualityIn1GVGP, slug


@testset "SPGainEqualityOutVGGP" begin
    @test SPGainEqualityOutVGGP <: SumProductRule{GainEquality}
    @test outboundType(SPGainEqualityOutVGGP) == Message{GaussianMeanPrecision}
    @test isApplicable(SPGainEqualityOutVGGP, [Nothing, Message{Gaussian}, Message{Gaussian}, Message{PointMass}])

    @test ruleSPGainEqualityOutVGGP(nothing, Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(2.0)), Message(Univariate, GaussianMeanPrecision, m=3.0, w=2.0), Message(Multivariate, PointMass, m=[1])) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0))
    @test ruleSPGainEqualityOutVGGP(nothing, Message(Multivariate, GaussianMeanPrecision, m=[0.0, 0.0], w=[1.0 0.0; 0.0 1.0]), Message(Univariate, GaussianMeanPrecision, m=3.0, w=2.0), Message(Multivariate, PointMass, m=[1.0, 1.0])) == Message(Multivariate, GaussianMeanPrecision, m=[1.2, 1.2], w=[3.0 2.0; 2.0 3.0])
    @test ruleSPGainEqualityOutVGGP(nothing, Message(Multivariate, GaussianMeanPrecision, m=[0.0, 0.0], w=[1.0 0.0; 0.0 1.0]), Message(Multivariate, GaussianMeanPrecision, m=[1.0, 1.0], w=[1.0 0.0; 0.0 1.0]), Message(MatrixVariate, PointMass, m=[1.0 1.0; 1.0 1.0])) == Message(Multivariate, GaussianMeanPrecision, m=[0.4, 0.4], w=[3.0 2.0; 2.0 3.0])
end

@testset "SPGainEqualityIn1GVGP" begin
    @test SPGainEqualityIn1GVGP <: SumProductRule{GainEquality}
    @test outboundType(SPGainEqualityIn1GVGP) == Message{GaussianMeanPrecision}
    @test isApplicable(SPGainEqualityIn1GVGP, [Message{Gaussian}, Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPGainEqualityIn1GVGP(Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(2.0)), nothing, Message(Univariate, GaussianMeanPrecision, m=3.0, w=2.0), Message(Multivariate, PointMass, m=[1])) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0))
end

@testset "slug" begin
    @test slug(GainEquality) == "gain_eq"
end

end # module
