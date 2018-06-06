module MultiplicationTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPMultiplicationOutVGP, SPMultiplicationOutVPG, SPMultiplicationIn1GVP, SPMultiplicationAGPV, SPMultiplicationOutVPP, SPMultiplicationIn1PVP, SPMultiplicationAPPV

@testset "Multiplication node construction through * syntax" begin
    g = FactorGraph()

    a = constant(1.0)
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = a*x
    @test isa(z, Variable)
    @test isa(g.nodes[:multiplication_1], Multiplication)

    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = 1.0*x
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
    @test ruleSPMultiplicationOutVGP(nothing, Message(Univariate, Gaussian, m=1.0, v=3.0), Message(Multivariate, PointMass, m=[2.0])) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutVGP(nothing, Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), Message(Univariate, PointMass, m=2.0)) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutVGP(nothing, Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutVPG" begin
    @test SPMultiplicationOutVPG <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPG) == Message{Gaussian}
    @test isApplicable(SPMultiplicationOutVPG, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPMultiplicationOutVPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, Gaussian, m=1.0, v=3.0)) == Message(Univariate, Gaussian, m=2.0, v=12.0)
    @test ruleSPMultiplicationOutVPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutVPG(nothing, Message(Multivariate, PointMass, m=[2.0]), Message(Univariate, Gaussian, m=1.0, v=3.0)) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutVPP" begin
    @test SPMultiplicationOutVPP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPP) == Message{PointMass}
    @test isApplicable(SPMultiplicationOutVPP, [Void, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPMultiplicationOutVPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPMultiplicationOutVPP(nothing, Message(Multivariate, PointMass, m=[2.0]), Message(MatrixVariate, PointMass, m=mat(1.0))) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPMultiplicationIn1GVP" begin
    @test SPMultiplicationIn1GVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1GVP) == Message{Gaussian}
    @test isApplicable(SPMultiplicationIn1GVP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1GVP(Message(Univariate, Gaussian, m=1.0, v=3.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=0.5, v=0.75)
    # @test ruleSPMultiplicationIn1GVP(Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), nothing, Message(Multivariate, PointMass, m=[2.0])) == Message(Univariate, Gaussian, m=0.5, v=0.75)
    @test ruleSPMultiplicationIn1GVP(Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[0.49999999999962497], v=mat(0.7499999999994372))
end

@testset "SPMultiplicationIn1PVP" begin
    @test SPMultiplicationIn1PVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1PVP) == Message{PointMass}
    @test isApplicable(SPMultiplicationIn1PVP, [Message{PointMass}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1PVP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, PointMass, m=0.5)
    @test ruleSPMultiplicationIn1PVP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, PointMass, m=[0.5])
end

@testset "SPMultiplicationAGPV" begin
    @test SPMultiplicationAGPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationAGPV) == Message{Gaussian}
    @test isApplicable(SPMultiplicationAGPV, [Message{Gaussian}, Message{PointMass}, Void]) 

    @test ruleSPMultiplicationAGPV(Message(Univariate, Gaussian, m=1.0, v=3.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, Gaussian, m=0.5, v=0.75)
    # @test ruleSPMultiplicationAGPV(Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), Message(Multivariate, PointMass, m=[2.0]), nothing) == Message(Univariate, Gaussian, m=0.5, v=0.75)
end

@testset "SPMultiplicationAPPV" begin
    @test SPMultiplicationAPPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationAPPV) == Message{PointMass}
    @test isApplicable(SPMultiplicationAPPV, [Message{PointMass}, Message{PointMass}, Void]) 

    @test ruleSPMultiplicationAPPV(Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, PointMass, m=0.5)
end

end # module