module MultiplicationTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPMultiplicationOutVGP, SPMultiplicationOutVPG, SPMultiplicationIn1GVP, SPMultiplicationIn2GPV

@testset "Multiplication node construction through * syntax" begin
    g = FactorGraph()

    a = constant(1.0)
    x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    z = a*x

    @test isa(z, Variable)
    @test isa(g.nodes[:multiplication_1], Multiplication)
end


#-------------
# Update rules
#-------------

@testset "SPMultiplicationOutVGP" begin
    @test SPMultiplicationOutVGP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVGP) == Message{Gaussian}
    @test isApplicable(SPMultiplicationOutVGP, [Void, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPMultiplicationOutVGP(nothing, Message(Univariate, Gaussian, m=1.0, v=3.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=2.0, v=12.0)
end

@testset "SPMultiplicationOutVPG" begin
    @test SPMultiplicationOutVPG <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPG) == Message{Gaussian}
    @test isApplicable(SPMultiplicationOutVPG, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPMultiplicationOutVPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, Gaussian, m=1.0, v=3.0)) == Message(Univariate, Gaussian, m=2.0, v=12.0)
end

@testset "SPMultiplicationIn1GVP" begin
    @test SPMultiplicationIn1GVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1GVP) == Message{Gaussian}
    @test isApplicable(SPMultiplicationIn1GVP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1GVP(Message(Univariate, Gaussian, m=1.0, v=3.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=0.5, v=0.75)
end

@testset "SPMultiplicationIn2GPV" begin
    @test SPMultiplicationIn2GPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn2GPV) == Message{Gaussian}
    @test isApplicable(SPMultiplicationIn2GPV, [Message{Gaussian}, Message{PointMass}, Void]) 

    @test ruleSPMultiplicationIn2GPV(Message(Univariate, Gaussian, m=1.0, v=3.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, Gaussian, m=0.5, v=0.75)
end

end # module