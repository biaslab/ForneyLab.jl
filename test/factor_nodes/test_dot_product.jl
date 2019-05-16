module DotProductTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPDotProductOutNPG, SPDotProductOutNGP, SPDotProductIn1GNP, SPDotProductIn2GPN


#-------------
# Update rules
#-------------

@testset "SPDotProductOutNPG" begin
    @test SPDotProductOutNPG <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutNPG) == Message{GaussianMeanVariance}
    @test isApplicable(SPDotProductOutNPG, [Nothing, Message{PointMass}, Message{Gaussian}])

    @test ruleSPDotProductOutNPG(nothing, Message(Multivariate, PointMass, m=[4.0]), Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(2.0))) == Message(Univariate, GaussianMeanVariance, m=12.0, v=32.0)
end

@testset "SPDotProductOutNGP" begin
    @test SPDotProductOutNGP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPDotProductOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPDotProductOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(2.0)), Message(Multivariate, PointMass, m=[4.0])) == Message(Univariate, GaussianMeanVariance, m=12.0, v=32.0)
end

@testset "SPDotProductIn1GNP" begin
    @test SPDotProductIn1GNP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn1GNP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPDotProductIn1GNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPDotProductIn1GNP(Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=2.0), nothing, Message(Multivariate, PointMass, m=[4.0])) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[12.0], w=mat(32.000000000001))
end

@testset "SPDotProductIn2GPN" begin
    @test SPDotProductIn2GPN <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn2GPN) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPDotProductIn2GPN, [Message{Gaussian}, Message{PointMass}, Nothing])

    @test ruleSPDotProductIn2GPN(Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=2.0), Message(Multivariate, PointMass, m=[4.0]), nothing) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[12.0], w=mat(32.000000000001))
end

end # module