module MultiplicationTest

using Test
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
    @test outboundType(SPMultiplicationOutVGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPMultiplicationOutVGP, [Nothing, Message{Gaussian}, Message{PointMass}]) 

    @test ruleSPMultiplicationOutVGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=12.0)
    @test ruleSPMultiplicationOutVGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0), Message(Multivariate, PointMass, m=[2.0])) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutVGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0)), Message(Univariate, PointMass, m=2.0)) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutVGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutVPG" begin
    @test SPMultiplicationOutVPG <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPG) == Message{GaussianMeanVariance}
    @test isApplicable(SPMultiplicationOutVPG, [Nothing, Message{PointMass}, Message{Gaussian}]) 

    @test ruleSPMultiplicationOutVPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=12.0)
    @test ruleSPMultiplicationOutVPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutVPG(nothing, Message(Multivariate, PointMass, m=[2.0]), Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutVPP" begin
    @test SPMultiplicationOutVPP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutVPP) == Message{PointMass}
    @test isApplicable(SPMultiplicationOutVPP, [Nothing, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPMultiplicationOutVPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPMultiplicationOutVPP(nothing, Message(Multivariate, PointMass, m=[2.0]), Message(MatrixVariate, PointMass, m=mat(1.0))) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPMultiplicationIn1GVP" begin
    @test SPMultiplicationIn1GVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1GVP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPMultiplicationIn1GVP, [Message{Gaussian}, Nothing, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1GVP(Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=3.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=12.0)
    @test ruleSPMultiplicationIn1GVP(Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(12.0 + tiny))
end

@testset "SPMultiplicationIn1PVP" begin
    @test SPMultiplicationIn1PVP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1PVP) == Message{PointMass}
    @test isApplicable(SPMultiplicationIn1PVP, [Message{PointMass}, Nothing, Message{PointMass}]) 

    @test ruleSPMultiplicationIn1PVP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, PointMass, m=0.5)
    @test ruleSPMultiplicationIn1PVP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, PointMass, m=[0.5])
end

@testset "SPMultiplicationAGPV" begin
    @test SPMultiplicationAGPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationAGPV) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPMultiplicationAGPV, [Message{Gaussian}, Message{PointMass}, Nothing]) 

    @test ruleSPMultiplicationAGPV(Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=3.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=12.0)
end

@testset "SPMultiplicationAPPV" begin
    @test SPMultiplicationAPPV <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationAPPV) == Message{PointMass}
    @test isApplicable(SPMultiplicationAPPV, [Message{PointMass}, Message{PointMass}, Nothing]) 

    @test ruleSPMultiplicationAPPV(Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, PointMass, m=0.5)
end

end # module