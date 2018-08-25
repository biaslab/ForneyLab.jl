module ExponentialTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPExponentialOutVG, SPExponentialOutVP, SPExponentialIn1LV, SPExponentialIn1PV

@testset "Exponential node construction through exp() syntax" begin
    g = FactorGraph()

    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = exp(x)

    @test isa(z, Variable)
    @test isa(g.nodes[:exponential_1], Exponential)
end


#-------------
# Update rules
#-------------

@testset "SPExponentialOutVG" begin
    @test SPExponentialOutVG <: SumProductRule{Exponential}
    @test outboundType(SPExponentialOutVG) == Message{LogNormal}
    @test isApplicable(SPExponentialOutVG, [Nothing, Message{Gaussian}])

    @test ruleSPExponentialOutVG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, LogNormal, m=1.0, s=2.0)
end

@testset "SPExponentialOutVP" begin
    @test SPExponentialOutVP <: SumProductRule{Exponential}
    @test outboundType(SPExponentialOutVP) == Message{PointMass}
    @test isApplicable(SPExponentialOutVP, [Nothing, Message{PointMass}])

    @test ruleSPExponentialOutVP(nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, PointMass, m=exp(2.0))
end

@testset "SPExponentialIn1LV" begin
    @test SPExponentialIn1LV <: SumProductRule{Exponential}
    @test outboundType(SPExponentialIn1LV) == Message{GaussianMeanVariance}
    @test isApplicable(SPExponentialIn1LV, [Message{LogNormal}, Nothing])

    @test ruleSPExponentialIn1LV(Message(Univariate, LogNormal, m=1.0, s=2.0), nothing) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
end

@testset "SPExponentialIn1PV" begin
    @test SPExponentialIn1PV <: SumProductRule{Exponential}
    @test outboundType(SPExponentialIn1PV) == Message{PointMass}
    @test isApplicable(SPExponentialIn1PV, [Message{PointMass}, Nothing])

    @test ruleSPExponentialIn1PV(Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, PointMass, m=log(2.0))
end

end # module