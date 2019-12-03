module ExponentialTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPExponentialOutNG, SPExponentialOutNP, SPExponentialIn1LN, SPExponentialIn1PN

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

@testset "SPExponentialOutNG" begin
    @test SPExponentialOutNG <: SumProductRule{Exponential}
    @test outboundType(SPExponentialOutNG) == Message{LogNormal}
    @test isApplicable(SPExponentialOutNG, [Nothing, Message{Gaussian}])

    @test ruleSPExponentialOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, LogNormal, m=1.0, s=2.0)
end

@testset "SPExponentialOutNP" begin
    @test SPExponentialOutNP <: SumProductRule{Exponential}
    @test outboundType(SPExponentialOutNP) == Message{PointMass}
    @test isApplicable(SPExponentialOutNP, [Nothing, Message{PointMass}])

    @test ruleSPExponentialOutNP(nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, PointMass, m=exp(2.0))
end

@testset "SPExponentialIn1LN" begin
    @test SPExponentialIn1LN <: SumProductRule{Exponential}
    @test outboundType(SPExponentialIn1LN) == Message{GaussianMeanVariance}
    @test isApplicable(SPExponentialIn1LN, [Message{LogNormal}, Nothing])

    @test ruleSPExponentialIn1LN(Message(Univariate, LogNormal, m=1.0, s=2.0), nothing) == Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)
end

@testset "SPExponentialIn1PN" begin
    @test SPExponentialIn1PN <: SumProductRule{Exponential}
    @test outboundType(SPExponentialIn1PN) == Message{PointMass}
    @test isApplicable(SPExponentialIn1PN, [Message{PointMass}, Nothing])

    @test ruleSPExponentialIn1PN(Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, PointMass, m=log(2.0))
end

end # module