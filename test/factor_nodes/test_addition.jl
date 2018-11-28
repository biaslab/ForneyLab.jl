module AdditionTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPAdditionOutVGG, SPAdditionOutVGP, SPAdditionOutVPG, SPAdditionIn1GVG, SPAdditionIn1PVG, SPAdditionIn2GGV, SPAdditionIn2PGV, SPAdditionIn1GVP, SPAdditionIn2GPV, SPAdditionOutVPP, SPAdditionIn1PVP, SPAdditionIn2PPV

@testset "Addition node construction through + syntax" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV y ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x + y
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)

    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x + 1.0
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
end

@testset "Addition node construction through - syntax" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV y ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x - y
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
    # This syntax works by changing the order interfaces are attached in. Below checks whether the pairings are succesful by
    # matching out and in2 interfaces of the corresponding nodes
    @test g.nodes[:addition_1].i[:out] == g.nodes[:gaussianmeanvariance_1].i[:out].partner
    @test g.nodes[:addition_1].i[:in2] == g.nodes[:gaussianmeanvariance_2].i[:out].partner

    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x - 1.0
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
end

#-------------
# Update rules
#-------------

@testset "SPAdditionOutVGG" begin
    @test SPAdditionOutVGG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVGG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionOutVGG, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(SPAdditionOutVGG, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPAdditionOutVGG, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleSPAdditionOutVGG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0)) == Message(Univariate, GaussianMeanVariance, m=4.0, v=6.0)
    @test ruleSPAdditionOutVGG(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0))) == Message(Multivariate, GaussianMeanVariance, m=[4.0], v=mat(6.0))
end

@testset "SPAdditionIn2GGV" begin
    @test SPAdditionIn2GGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GGV) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn2GGV, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleSPAdditionIn2GGV(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.0, v=6.0)
    @test ruleSPAdditionIn2GGV(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(6.0))
end

@testset "SPAdditionIn1GVG" begin
    @test SPAdditionIn1GVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GVG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn1GVG, [Message{Gaussian}, Nothing, Message{Gaussian}])

    @test ruleSPAdditionIn1GVG(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=6.0)
    @test ruleSPAdditionIn1GVG(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(6.0))
end

@testset "SPAdditionOutVGP" begin
    @test SPAdditionOutVGP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionOutVGP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPAdditionOutVGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
    @test ruleSPAdditionOutVGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), Message(Multivariate, PointMass, m=[3.0])) == Message(Multivariate, GaussianMeanVariance, m=[4.0], v=mat(2.0))
end

@testset "SPAdditionOutVPG" begin
    @test SPAdditionOutVPG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVPG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionOutVPG, [Nothing, Message{PointMass}, Message{Gaussian}])

    @test ruleSPAdditionOutVPG(nothing, Message(Univariate, PointMass, m=3.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
    @test ruleSPAdditionOutVPG(nothing, Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[4.0], v=mat(2.0))
end

@testset "SPAdditionIn1PVG" begin
    @test SPAdditionIn1PVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PVG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn1PVG, [Message{PointMass}, Nothing, Message{Gaussian}])

    @test ruleSPAdditionIn1PVG(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0)
    @test ruleSPAdditionIn1PVG(Message(Multivariate, PointMass, m=[3.0]), nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(2.0))
end

@testset "SPAdditionIn2PGV" begin
    @test SPAdditionIn2PGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PGV) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn2PGV, [Message{PointMass}, Message{Gaussian}, Nothing])

    @test ruleSPAdditionIn2PGV(Message(Univariate, PointMass, m=3.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0)
    @test ruleSPAdditionIn2PGV(Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(2.0))
end

@testset "SPAdditionIn1GVP" begin
    @test SPAdditionIn1GVP <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GVP) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn1GVP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPAdditionIn1GVP(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=4.0)
    @test ruleSPAdditionIn1GVP(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))
end

@testset "SPAdditionIn2GPV" begin
    @test SPAdditionIn2GPV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GPV) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn2GPV, [Message{Gaussian}, Message{PointMass}, Nothing])

    @test ruleSPAdditionIn2GPV(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.0, v=4.0)
    @test ruleSPAdditionIn2GPV(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Message(Multivariate, PointMass, m=[1.0]), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))
end

@testset "SPAdditionOutVPP" begin
    @test SPAdditionOutVPP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVPP) == Message{PointMass}
    @test isApplicable(SPAdditionOutVPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPAdditionOutVPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, PointMass, m=4.0)
    @test ruleSPAdditionOutVPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(Multivariate, PointMass, m=[3.0])) == Message(Multivariate, PointMass, m=[4.0])
end

@testset "SPAdditionIn2PPV" begin
    @test SPAdditionIn2PPV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PPV) == Message{PointMass}
    @test isApplicable(SPAdditionIn2PPV, [Message{PointMass}, Message{PointMass}, Nothing])

    @test ruleSPAdditionIn2PPV(Message(Univariate, PointMass, m=3.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPAdditionIn2PPV(Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, PointMass, m=[1.0]), nothing) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPAdditionIn1PVP" begin
    @test SPAdditionIn1PVP <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PVP) == Message{PointMass}
    @test isApplicable(SPAdditionIn1PVP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPAdditionIn1PVP(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPAdditionIn1PVP(Message(Multivariate, PointMass, m=[3.0]), nothing, Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, PointMass, m=[2.0])
end

end # module
