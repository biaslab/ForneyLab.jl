module DotProductTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPDotProductOutVPG, SPDotProductOutVGP, SPDotProductIn1GVP, SPDotProductIn2GPV


#-------------
# Update rules
#-------------

@testset "SPDotProductOutVPG" begin
    @test SPDotProductOutVPG <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutVPG) == Message{GaussianMeanVariance}
    @test isApplicable(SPDotProductOutVPG, [Nothing, Message{PointMass}, Message{Gaussian}])

    @test ruleSPDotProductOutVPG(nothing, Message(Multivariate, PointMass, m=[4.0]), Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(2.0))) == Message(Univariate, GaussianMeanVariance, m=12.0, v=32.0)
end

@testset "SPDotProductOutVGP" begin
    @test SPDotProductOutVGP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutVGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPDotProductOutVGP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPDotProductOutVGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(2.0)), Message(Multivariate, PointMass, m=[4.0])) == Message(Univariate, GaussianMeanVariance, m=12.0, v=32.0)
end

@testset "SPDotProductIn1GVP" begin
    @test SPDotProductIn1GVP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn1GVP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPDotProductIn1GVP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPDotProductIn1GVP(Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=2.0), nothing, Message(Multivariate, PointMass, m=[4.0])) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[12.0], w=mat(32.000000000001))
end

@testset "SPDotProductIn2GPV" begin
    @test SPDotProductIn2GPV <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn2GPV) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPDotProductIn2GPV, [Message{Gaussian}, Message{PointMass}, Nothing])

    @test ruleSPDotProductIn2GPV(Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=2.0), Message(Multivariate, PointMass, m=[4.0]), nothing) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[12.0], w=mat(32.000000000001))
end

end # module