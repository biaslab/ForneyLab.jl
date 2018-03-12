module DotProductTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPDotProductOutVPG, SPDotProductOutVGP, SPDotProductIn1GVP, SPDotProductIn2GPV


#-------------
# Update rules
#-------------

@testset "SPDotProductOutVPG" begin
    @test SPDotProductOutVPG <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutVPG) == Message{Gaussian}
    @test isApplicable(SPDotProductOutVPG, [Void, Message{PointMass}, Message{Gaussian}])

    @test ruleSPDotProductOutVPG(nothing, Message(Multivariate, PointMass, m=[4.0]), Message(Multivariate, Gaussian, m=[3.0], v=mat(2.0))) == Message(Univariate, Gaussian, m=12.0, v=32.0)
end

@testset "SPDotProductOutVGP" begin
    @test SPDotProductOutVGP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductOutVGP) == Message{Gaussian}
    @test isApplicable(SPDotProductOutVGP, [Void, Message{Gaussian}, Message{PointMass}])

    @test ruleSPDotProductOutVGP(nothing, Message(Multivariate, Gaussian, m=[3.0], v=mat(2.0)), Message(Multivariate, PointMass, m=[4.0])) == Message(Univariate, Gaussian, m=12.0, v=32.0)
end

@testset "SPDotProductIn1GVP" begin
    @test SPDotProductIn1GVP <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn1GVP) == Message{Gaussian}
    @test isApplicable(SPDotProductIn1GVP, [Message{Gaussian}, Void, Message{PointMass}])

    @test ruleSPDotProductIn1GVP(Message(Univariate, Gaussian, xi=3.0, w=2.0), nothing, Message(Multivariate, PointMass, m=[4.0])) == Message(Multivariate, Gaussian, xi=[12.0], w=mat(32.0))
end

@testset "SPDotProductIn2GPV" begin
    @test SPDotProductIn2GPV <: SumProductRule{DotProduct}
    @test outboundType(SPDotProductIn2GPV) == Message{Gaussian}
    @test isApplicable(SPDotProductIn2GPV, [Message{Gaussian}, Message{PointMass}, Void])

    @test ruleSPDotProductIn2GPV(Message(Univariate, Gaussian, xi=3.0, w=2.0), Message(Multivariate, PointMass, m=[4.0]), nothing) == Message(Multivariate, Gaussian, xi=[12.0], w=mat(32.0))
end

end # module