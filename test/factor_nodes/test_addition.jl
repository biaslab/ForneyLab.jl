module AdditionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPAdditionOutVGG, SPAdditionOutVGP, SPAdditionOutVPG, SPAdditionIn1GVG, SPAdditionIn1PVG, SPAdditionIn2GGV, SPAdditionIn2PGV

@testset "Addition node construction through + syntax" begin
    g = FactorGraph()

    x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    y ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    z = x + y

    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
end


#-------------
# Update rules
#-------------

@testset "SPAdditionOutVGG" begin
    @test SPAdditionOutVGG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVGG) == Message{Gaussian}
    @test isApplicable(SPAdditionOutVGG, [Void, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(SPAdditionOutVGG, [Void, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPAdditionOutVGG, [Message{Gaussian}, Message{Gaussian}, Void])

    @test ruleSPAdditionOutVGG(nothing, Message(Univariate, Gaussian, m=1.0, v=2.0), Message(Univariate, Gaussian, m=3.0, v=4.0)) == Message(Univariate, Gaussian, m=4.0, v=6.0)
end

@testset "SPAdditionIn2GGV" begin
    @test SPAdditionIn2GGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GGV) == Message{Gaussian}
    @test isApplicable(SPAdditionIn2GGV, [Message{Gaussian}, Message{Gaussian}, Void])

    @test ruleSPAdditionIn2GGV(Message(Univariate, Gaussian, m=3.0, v=4.0), Message(Univariate, Gaussian, m=1.0, v=2.0), nothing) == Message(Univariate, Gaussian, m=2.0, v=6.0)
end

@testset "SPAdditionIn1GVG" begin
    @test SPAdditionIn1GVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GVG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1GVG, [Message{Gaussian}, Void, Message{Gaussian}])

    @test ruleSPAdditionIn1GVG(Message(Univariate, Gaussian, m=3.0, v=4.0), nothing, Message(Univariate, Gaussian, m=1.0, v=2.0)) == Message(Univariate, Gaussian, m=2.0, v=6.0)
end

@testset "SPAdditionOutVGP" begin
    @test SPAdditionOutVGP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVGP) == Message{Gaussian}
    @test isApplicable(SPAdditionOutVGP, [Void, Message{Gaussian}, Message{PointMass}])

    @test ruleSPAdditionOutVGP(nothing, Message(Univariate, Gaussian, m=1.0, v=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian, m=4.0, v=2.0)
end

@testset "SPAdditionOutVPG" begin
    @test SPAdditionOutVPG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVPG) == Message{Gaussian}
    @test isApplicable(SPAdditionOutVPG, [Void, Message{PointMass}, Message{Gaussian}])

    @test ruleSPAdditionOutVPG(nothing, Message(Univariate, PointMass, m=3.0), Message(Univariate, Gaussian, m=1.0, v=2.0)) == Message(Univariate, Gaussian, m=4.0, v=2.0)
end

@testset "SPAdditionIn1PVG" begin
    @test SPAdditionIn1PVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PVG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1PVG, [Message{PointMass}, Void, Message{Gaussian}])

    @test ruleSPAdditionIn1PVG(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, Gaussian, m=1.0, v=2.0)) == Message(Univariate, Gaussian, m=2.0, v=2.0)
end

@testset "SPAdditionIn2PGV" begin
    @test SPAdditionIn2PGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PGV) == Message{Gaussian}
    @test isApplicable(SPAdditionIn2PGV, [Message{PointMass}, Message{Gaussian}, Void])

    @test ruleSPAdditionIn2PGV(Message(Univariate, PointMass, m=3.0), Message(Univariate, Gaussian, m=1.0, v=2.0), nothing) == Message(Univariate, Gaussian, m=2.0, v=2.0)
end

end # module