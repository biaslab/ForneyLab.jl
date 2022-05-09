module DotProductTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPDotProductOutNPG, SPDotProductOutNGP, SPDotProductIn1GNP, SPDotProductIn2GPN


#-------------
# Update rules
#-------------

@testset "SPDotProductOutNPG" begin
    @test SPDotProductOutNPG <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutNPG) == Message{Gaussian{Moments}}
    @test isApplicable(SPDotProductOutNPG, [Nothing, Message{PointMass}, Message{Gaussian}])

    @test ruleSPDotProductOutNPG(nothing, Message(Multivariate, PointMass, m=[4.0]), Message(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(2.0))) == Message(Univariate, Gaussian{Moments}, m=12.0, v=32.0)
end

@testset "SPDotProductOutNGP" begin
    @test SPDotProductOutNGP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutNGP) == Message{Gaussian{Moments}}
    @test isApplicable(SPDotProductOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPDotProductOutNGP(nothing, Message(Multivariate, Gaussian{Moments}, m=[3.0], v=mat(2.0)), Message(Multivariate, PointMass, m=[4.0])) == Message(Univariate, Gaussian{Moments}, m=12.0, v=32.0)
end

@testset "SPDotProductIn1GNP" begin
    @test SPDotProductIn1GNP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn1GNP) == Message{Gaussian{Canonical}}
    @test isApplicable(SPDotProductIn1GNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPDotProductIn1GNP(Message(Univariate, Gaussian{Canonical}, xi=3.0, w=2.0), nothing, Message(Multivariate, PointMass, m=[4.0])) == Message(Multivariate, Gaussian{Canonical}, xi=[12.0], w=mat(32.000000000001))
end

@testset "SPDotProductIn2GPN" begin
    @test SPDotProductIn2GPN <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn2GPN) == Message{Gaussian{Canonical}}
    @test isApplicable(SPDotProductIn2GPN, [Message{Gaussian}, Message{PointMass}, Nothing])

    @test ruleSPDotProductIn2GPN(Message(Univariate, Gaussian{Canonical}, xi=3.0, w=2.0), Message(Multivariate, PointMass, m=[4.0]), nothing) == Message(Multivariate, Gaussian{Canonical}, xi=[12.0], w=mat(32.000000000001))
end

end # module