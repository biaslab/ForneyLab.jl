module MultiplicationTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPMultiplicationOutGP, SPMultiplicationInGP

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

@testset "SPMultiplicationOutGP" begin
    @test SPMultiplicationOutGP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutGP) == Message{Gaussian}
    @test isApplicable(SPMultiplicationOutGP, [Void, Message{Gaussian}, Message{PointMass}]) 
    @test !isApplicable(SPMultiplicationOutGP, [Void, Message{PointMass}, Message{Gaussian}]) 
    @test !isApplicable(SPMultiplicationOutGP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationOutGP(nothing, Message(Univariate(Gaussian, m=1.0, v=3.0)), Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=2.0, v=12.0))
end

@testset "SPMultiplicationInGP" begin
    @test SPMultiplicationInGP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationInGP) == Message{Gaussian}
    @test !isApplicable(SPMultiplicationInGP, [Void, Message{Gaussian}, Message{PointMass}]) 
    @test !isApplicable(SPMultiplicationInGP, [Message{PointMass}, Void, Message{Gaussian}]) 
    @test isApplicable(SPMultiplicationInGP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationInGP(Message(Univariate(Gaussian, m=1.0, v=3.0)), nothing, Message(Univariate(PointMass, m=2.0))) == Message(Univariate(Gaussian, m=0.5, v=0.75))
end

end # module