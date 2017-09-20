module MultiplicationTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPMultiplicationGPV, SPMultiplicationVPG

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

@testset "SPMultiplicationGPV" begin
    @test SPMultiplicationGPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationGPV) == Message{Gaussian}
    @test isApplicable(SPMultiplicationGPV, [Message{Gaussian}, Message{PointMass}, Void]) 
    @test !isApplicable(SPMultiplicationGPV, [Message{PointMass}, Message{Gaussian}, Void]) 
    @test !isApplicable(SPMultiplicationGPV, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPMultiplicationGPV(Message(Gaussian, m=1.0, v=3.0), Message(PointMass, m=2.0), nothing) == Message(Gaussian, m=2.0, v=12.0)
end

@testset "SPMultiplicationVPG" begin
    @test SPMultiplicationVPG <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationVPG) == Message{Gaussian}
    @test !isApplicable(SPMultiplicationVPG, [Message{Gaussian}, Message{PointMass}, Void]) 
    @test !isApplicable(SPMultiplicationVPG, [Void, Message{Gaussian}, Message{PointMass}]) 
    @test isApplicable(SPMultiplicationVPG, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPMultiplicationVPG(nothing, Message(PointMass, m=2.0), Message(Gaussian, m=1.0, v=3.0)) == Message(Gaussian, m=0.5, v=0.75)
end

end # module