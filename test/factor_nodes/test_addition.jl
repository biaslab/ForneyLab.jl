module AdditionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPAdditionOutVGG, SPAdditionOutVGP, SPAdditionOutVPG, SPAdditionIn1GVG, SPAdditionIn1PVG, SPAdditionIn2GGV, SPAdditionIn2PGV, SPAdditionIn1GVP, SPAdditionIn2GPV, SPAdditionOutVPP, SPAdditionIn1PVP, SPAdditionIn2PPV

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
    @test ruleSPAdditionOutVGG(nothing, Message(Multivariate, Gaussian, m=[1.0], v=[2.0].'), Message(Multivariate, Gaussian, m=[3.0], v=[4.0].')) == Message(Multivariate, Gaussian, m=[4.0], v=[6.0].')
end

@testset "SPAdditionIn2GGV" begin
    @test SPAdditionIn2GGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GGV) == Message{Gaussian}
    @test isApplicable(SPAdditionIn2GGV, [Message{Gaussian}, Message{Gaussian}, Void])

    @test ruleSPAdditionIn2GGV(Message(Univariate, Gaussian, m=3.0, v=4.0), Message(Univariate, Gaussian, m=1.0, v=2.0), nothing) == Message(Univariate, Gaussian, m=2.0, v=6.0)
    @test ruleSPAdditionIn2GGV(Message(Multivariate, Gaussian, m=[3.0], v=[4.0].'), Message(Multivariate, Gaussian, m=[1.0], v=[2.0].'), nothing) == Message(Multivariate, Gaussian, m=[2.0], v=[6.0].')
end

@testset "SPAdditionIn1GVG" begin
    @test SPAdditionIn1GVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GVG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1GVG, [Message{Gaussian}, Void, Message{Gaussian}])

    @test ruleSPAdditionIn1GVG(Message(Univariate, Gaussian, m=3.0, v=4.0), nothing, Message(Univariate, Gaussian, m=1.0, v=2.0)) == Message(Univariate, Gaussian, m=2.0, v=6.0)
    @test ruleSPAdditionIn1GVG(Message(Multivariate, Gaussian, m=[3.0], v=[4.0].'), nothing, Message(Multivariate, Gaussian, m=[1.0], v=[2.0].')) == Message(Multivariate, Gaussian, m=[2.0], v=[6.0].')
end

@testset "SPAdditionOutVGP" begin
    @test SPAdditionOutVGP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVGP) == Message{Gaussian}
    @test isApplicable(SPAdditionOutVGP, [Void, Message{Gaussian}, Message{PointMass}])

    @test ruleSPAdditionOutVGP(nothing, Message(Univariate, Gaussian, m=1.0, v=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, Gaussian, m=4.0, v=2.0)
    @test ruleSPAdditionOutVGP(nothing, Message(Multivariate, Gaussian, m=[1.0], v=[2.0].'), Message(Multivariate, PointMass, m=[3.0])) == Message(Multivariate, Gaussian, m=[4.0], v=[2.0].')
end

@testset "SPAdditionOutVPG" begin
    @test SPAdditionOutVPG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVPG) == Message{Gaussian}
    @test isApplicable(SPAdditionOutVPG, [Void, Message{PointMass}, Message{Gaussian}])

    @test ruleSPAdditionOutVPG(nothing, Message(Univariate, PointMass, m=3.0), Message(Univariate, Gaussian, m=1.0, v=2.0)) == Message(Univariate, Gaussian, m=4.0, v=2.0)
    @test ruleSPAdditionOutVPG(nothing, Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, Gaussian, m=[1.0], v=[2.0].')) == Message(Multivariate, Gaussian, m=[4.0], v=[2.0].')
end

@testset "SPAdditionIn1PVG" begin
    @test SPAdditionIn1PVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PVG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1PVG, [Message{PointMass}, Void, Message{Gaussian}])

    @test ruleSPAdditionIn1PVG(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, Gaussian, m=1.0, v=2.0)) == Message(Univariate, Gaussian, m=2.0, v=2.0)
    @test ruleSPAdditionIn1PVG(Message(Multivariate, PointMass, m=[3.0]), nothing, Message(Multivariate, Gaussian, m=[1.0], v=[2.0].')) == Message(Multivariate, Gaussian, m=[2.0], v=[2.0].')
end

@testset "SPAdditionIn2PGV" begin
    @test SPAdditionIn2PGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PGV) == Message{Gaussian}
    @test isApplicable(SPAdditionIn2PGV, [Message{PointMass}, Message{Gaussian}, Void])

    @test ruleSPAdditionIn2PGV(Message(Univariate, PointMass, m=3.0), Message(Univariate, Gaussian, m=1.0, v=2.0), nothing) == Message(Univariate, Gaussian, m=2.0, v=2.0)
    @test ruleSPAdditionIn2PGV(Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, Gaussian, m=[1.0], v=[2.0].'), nothing) == Message(Multivariate, Gaussian, m=[2.0], v=[2.0].')
end

@testset "SPAdditionIn1GVP" begin
    @test SPAdditionIn1GVP <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GVP) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1GVP, [Message{Gaussian}, Void, Message{PointMass}])

    @test ruleSPAdditionIn1GVP(Message(Univariate, Gaussian, m=3.0, v=4.0), nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, Gaussian, m=2.0, v=4.0)
    @test ruleSPAdditionIn1GVP(Message(Multivariate, Gaussian, m=[3.0], v=[4.0].'), nothing, Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, Gaussian, m=[2.0], v=[4.0].')
end

@testset "SPAdditionIn2GPV" begin
    @test SPAdditionIn2GPV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GPV) == Message{Gaussian}
    @test isApplicable(SPAdditionIn2GPV, [Message{Gaussian}, Message{PointMass}, Void])

    @test ruleSPAdditionIn2GPV(Message(Univariate, Gaussian, m=3.0, v=4.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, Gaussian, m=2.0, v=4.0)
    @test ruleSPAdditionIn2GPV(Message(Multivariate, Gaussian, m=[3.0], v=[4.0].'), Message(Multivariate, PointMass, m=[1.0]), nothing) == Message(Multivariate, Gaussian, m=[2.0], v=[4.0].')
end

@testset "SPAdditionOutVPP" begin
    @test SPAdditionOutVPP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutVPP) == Message{PointMass}
    @test isApplicable(SPAdditionOutVPP, [Void, Message{PointMass}, Message{PointMass}])

    @test ruleSPAdditionOutVPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, PointMass, m=4.0)
    @test ruleSPAdditionOutVPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(Multivariate, PointMass, m=[3.0])) == Message(Multivariate, PointMass, m=[4.0])
end

@testset "SPAdditionIn2PPV" begin
    @test SPAdditionIn2PPV <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PPV) == Message{PointMass}
    @test isApplicable(SPAdditionIn2PPV, [Message{PointMass}, Message{PointMass}, Void])

    @test ruleSPAdditionIn2PPV(Message(Univariate, PointMass, m=3.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPAdditionIn2PPV(Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, PointMass, m=[1.0]), nothing) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPAdditionIn1PVP" begin
    @test SPAdditionIn1PVP <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PVP) == Message{PointMass}
    @test isApplicable(SPAdditionIn1PVP, [Message{PointMass}, Void, Message{PointMass}])

    @test ruleSPAdditionIn1PVP(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPAdditionIn1PVP(Message(Multivariate, PointMass, m=[3.0]), nothing, Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, PointMass, m=[2.0])
end

end # module