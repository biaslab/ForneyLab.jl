module AdditionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPAdditionOutGG, SPAdditionOutGP, SPAdditionOutPG, SPAdditionIn1GG, SPAdditionIn1PG, SPAdditionIn2GG

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

@testset "SPAdditionOutGG" begin
    @test SPAdditionOutGG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutGG) == Message{Gaussian}
    @test isApplicable(SPAdditionOutGG, [Void, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(SPAdditionOutGG, [Void, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPAdditionOutGG, [Message{Gaussian}, Message{Gaussian}, Void])

    @test ruleSPAdditionOutGG(nothing, Message(Univariate(Gaussian, m=1.0, v=2.0)), Message(Univariate(Gaussian, m=3.0, v=4.0))) == Message(Univariate(Gaussian, m=4.0, v=6.0))
end

@testset "SPAdditionIn2GG" begin
    @test SPAdditionIn2GG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn2GG, [Message{Gaussian}, Message{Gaussian}, Void])

    @test ruleSPAdditionIn2GG(Message(Univariate(Gaussian, m=3.0, v=4.0)), Message(Univariate(Gaussian, m=1.0, v=2.0)), nothing) == Message(Univariate(Gaussian, m=2.0, v=6.0))
end

@testset "SPAdditionIn1GG" begin
    @test SPAdditionIn1GG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1GG, [Message{Gaussian}, Void, Message{Gaussian}])

    @test ruleSPAdditionIn1GG(Message(Univariate(Gaussian, m=3.0, v=4.0)), nothing, Message(Univariate(Gaussian, m=1.0, v=2.0))) == Message(Univariate(Gaussian, m=2.0, v=6.0))
end

@testset "SPAdditionOutGP" begin
    @test SPAdditionOutGP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutGP) == Message{Gaussian}
    @test isApplicable(SPAdditionOutGP, [Void, Message{Gaussian}, Message{PointMass}])

    @test ruleSPAdditionOutGP(nothing, Message(Univariate(Gaussian, m=1.0, v=2.0)), Message(Univariate(PointMass, m=3.0))) == Message(Univariate(Gaussian, m=4.0, v=2.0))
end

@testset "SPAdditionOutPG" begin
    @test SPAdditionOutPG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutPG) == Message{Gaussian}
    @test isApplicable(SPAdditionOutPG, [Void, Message{PointMass}, Message{Gaussian}])

    @test ruleSPAdditionOutPG(nothing, Message(Univariate(PointMass, m=3.0)), Message(Univariate(Gaussian, m=1.0, v=2.0))) == Message(Univariate(Gaussian, m=4.0, v=2.0))
end
@testset "SPAdditionIn1PG" begin
    @test SPAdditionIn1PG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PG) == Message{Gaussian}
    @test isApplicable(SPAdditionIn1PG, [Message{PointMass}, Void, Message{Gaussian}])

    @test ruleSPAdditionIn1PG(Message(Univariate(PointMass, m=3.0)), nothing, Message(Univariate(Gaussian, m=1.0, v=2.0))) == Message(Univariate(Gaussian, m=2.0, v=2.0))
end

end # module