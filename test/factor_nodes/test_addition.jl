module AdditionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPAdditionGGV, SPAdditionGVG, SPAdditionVGG, SPAdditionVGP, SPAdditionGPV, SPAdditionPGV

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

@testset "SPAdditionGGV" begin
    @test SPAdditionGGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionGGV) == Message{Gaussian}
    @test isApplicable(SPAdditionGGV, [Message{Gaussian}, Message{Gaussian}, Void]) 
    @test !isApplicable(SPAdditionGGV, [Message{PointMass}, Message{PointMass}, Void]) 
    @test !isApplicable(SPAdditionGGV, [Void, Message{Gaussian}, Message{Gaussian}]) 

    @test ruleSPAdditionGGV(Message(Gaussian, m=1.0, v=2.0), Message(Gaussian, m=3.0, v=4.0), nothing) == Message(Gaussian, m=4.0, v=6.0)
end

@testset "SPAdditionGVG" begin
    @test SPAdditionGVG <: SumProductRule{Addition}
    @test outboundType(SPAdditionGVG) == Message{Gaussian}
    @test isApplicable(SPAdditionGVG, [Message{Gaussian}, Void, Message{Gaussian}]) 

    @test ruleSPAdditionGVG(Message(Gaussian, m=1.0, v=2.0), nothing, Message(Gaussian, m=3.0, v=4.0)) == Message(Gaussian, m=2.0, v=6.0)
end

@testset "SPAdditionVGG" begin
    @test SPAdditionVGG <: SumProductRule{Addition}
    @test outboundType(SPAdditionVGG) == Message{Gaussian}
    @test isApplicable(SPAdditionVGG, [Void, Message{Gaussian}, Message{Gaussian}]) 

    @test ruleSPAdditionVGG(nothing, Message(Gaussian, m=1.0, v=2.0), Message(Gaussian, m=3.0, v=4.0)) == Message(Gaussian, m=2.0, v=6.0)
end

@testset "SPAdditionGPV" begin
    @test SPAdditionGPV <: SumProductRule{Addition}
    @test outboundType(SPAdditionGPV) == Message{Gaussian}
    @test isApplicable(SPAdditionGPV, [Message{Gaussian}, Message{PointMass}, Void]) 

    @test ruleSPAdditionGPV(Message(Gaussian, m=1.0, v=2.0), Message(PointMass, m=3.0), nothing) == Message(Gaussian, m=4.0, v=2.0)
end

@testset "SPAdditionPGV" begin
    @test SPAdditionPGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionPGV) == Message{Gaussian}
    @test isApplicable(SPAdditionPGV, [Message{PointMass}, Message{Gaussian}, Void]) 

    @test ruleSPAdditionPGV(Message(PointMass, m=3.0), Message(Gaussian, m=1.0, v=2.0), nothing) == Message(Gaussian, m=4.0, v=2.0)
end
@testset "SPAdditionVGP" begin
    @test SPAdditionVGP <: SumProductRule{Addition}
    @test outboundType(SPAdditionVGP) == Message{Gaussian}
    @test isApplicable(SPAdditionVGP, [Void, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPAdditionVGP(nothing, Message(Gaussian, m=1.0, v=2.0), Message(PointMass, m=3.0)) == Message(Gaussian, m=2.0, v=2.0)
end

end # module