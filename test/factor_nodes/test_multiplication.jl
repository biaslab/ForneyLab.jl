module MultiplicationTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPMultiplicationOutVGP, SPMultiplicationInGVP

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
    @test !isApplicable(SPMultiplicationOutVGP, [Void, Message{PointMass}, Message{Gaussian}]) 
    @test !isApplicable(SPMultiplicationOutVGP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationOutVGP(nothing, Message(Univariate(Gaussian, m=1.0, v=3.0)), Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=2.0, v=12.0))
end

@testset "SPMultiplicationInGVP" begin
    @test SPMultiplicationInGVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationInGVP) == Message{Gaussian}
    @test !isApplicable(SPMultiplicationInGVP, [Void, Message{Gaussian}, Message{PointMass}]) 
    @test !isApplicable(SPMultiplicationInGVP, [Message{PointMass}, Void, Message{Gaussian}]) 
    @test isApplicable(SPMultiplicationInGVP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationInGVP(Message(Univariate(Gaussian, m=1.0, v=3.0)), nothing, Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=0.5, v=0.75))
end

end # module