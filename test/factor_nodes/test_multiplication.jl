module MultiplicationTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPMultiplicationOutVGP, SPMultiplicationOutVPG, SPMultiplicationIn1GVP, SPMultiplicationIn2GPV, SPMultiplicationOutVPP, SPMultiplicationIn1PVP, SPMultiplicationIn2PPV

@testset "Multiplication node construction through * syntax" begin
    g = FactorGraph()

    a = constant(1.0)
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = a*x

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
    @test ruleSPMultiplicationOutVGP(nothing, Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutVPG" begin
    @test SPMultiplicationOutVPG <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPG) == Message{Gaussian}
    @test isApplicable(SPMultiplicationOutVPG, [Void, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPMultiplicationOutVPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, Gaussian, m=1.0, v=3.0)) == Message(Univariate, Gaussian, m=2.0, v=12.0)
    @test ruleSPMultiplicationOutVPG(nothing, Message(MatrixVariate, PointMass, m=mat(2.0)), Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0))) == Message(Multivariate, Gaussian, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutVPP" begin
    @test SPMultiplicationOutVPP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPP) == Message{PointMass}
    @test isApplicable(SPMultiplicationOutVPP, [Void, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPMultiplicationOutVPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPMultiplicationOutVPP(nothing, Message(MatrixVariate, PointMass, m=mat(2.0)), Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, PointMass, m=[2.0])
    @test ruleSPMultiplicationOutVPP(nothing, Message(Multivariate, PointMass, m=[2.0]), Message(MatrixVariate, PointMass, m=mat(1.0))) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPMultiplicationIn1GVP" begin
    @test SPMultiplicationIn1GVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1GVP) == Message{Gaussian}
    @test isApplicable(SPMultiplicationIn1GVP, [Message{Gaussian}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1GVP(Message(Univariate, Gaussian, m=1.0, v=3.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, Gaussian, m=0.5, v=0.75)
    @test ruleSPMultiplicationIn1GVP(Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, Gaussian, m=[0.5], v=mat(0.75))
end

@testset "SPMultiplicationIn1PVP" begin
    @test SPMultiplicationIn1PVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1PVP) == Message{PointMass}
    @test isApplicable(SPMultiplicationIn1PVP, [Message{PointMass}, Void, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1PVP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, PointMass, m=0.5)
    @test ruleSPMultiplicationIn1PVP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, PointMass, m=[0.5])
end

@testset "SPMultiplicationIn2GPV" begin
    @test SPMultiplicationIn2GPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn2GPV) == Message{Gaussian}
    @test isApplicable(SPMultiplicationIn2GPV, [Message{Gaussian}, Message{PointMass}, Void]) 

    @test ruleSPMultiplicationIn2GPV(Message(Univariate, Gaussian, m=1.0, v=3.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, Gaussian, m=0.5, v=0.75)
    @test ruleSPMultiplicationIn2GPV(Message(Multivariate, Gaussian, m=[1.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0)), nothing) == Message(Multivariate, Gaussian, m=[0.5], v=mat(0.75))
end

@testset "SPMultiplicationIn2PPV" begin
    @test SPMultiplicationIn2PPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn2PPV) == Message{PointMass}
    @test isApplicable(SPMultiplicationIn2PPV, [Message{PointMass}, Message{PointMass}, Void]) 

    @test ruleSPMultiplicationIn2PPV(Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, PointMass, m=0.5)
    @test ruleSPMultiplicationIn2PPV(Message(Multivariate, PointMass, m=[1.0]), Message(MatrixVariate, PointMass, m=mat(2.0)), nothing) == Message(Multivariate, PointMass, m=[0.5])
end

end # module