module MultiplicationTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPMultiplicationOutNGP, SPMultiplicationOutNPG, SPMultiplicationIn1GNP, SPMultiplicationAGPN, SPMultiplicationOutNPP, SPMultiplicationIn1PNP, SPMultiplicationAPPN

@testset "Multiplication node construction through * syntax" begin
    g = FactorGraph()
    a = constant(1.0)
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = a*x
    @test isa(z, Variable)
    @test isa(g.nodes[:multiplication_1], Multiplication)

    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = 1.0*x
    @test isa(z, Variable)
    @test isa(g.nodes[:multiplication_1], Multiplication)

    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x*1.0
    @test isa(z, Variable)
    @test isa(g.nodes[:multiplication_1], Multiplication)
end


#-------------
# Update rules
#-------------

@testset "SPMultiplicationOutNGP" begin
    @test SPMultiplicationOutNGP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPMultiplicationOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPMultiplicationOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0), Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=12.0)
    @test ruleSPMultiplicationOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0), Message(Multivariate, PointMass, m=[2.0])) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
    @test ruleSPMultiplicationOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(3.0)), Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(12.0))
end

@testset "SPMultiplicationOutNPG" begin
    @test SPMultiplicationOutNPG <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutNPG) == Message{GaussianMeanVariance}
    @test isApplicable(SPMultiplicationOutNPG, [Nothing, Message{PointMass}, Message{Gaussian}])

    @test ruleSPMultiplicationOutNPG(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=12.0)
end

@testset "SPMultiplicationOutNPP" begin
    @test SPMultiplicationOutNPP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationOutNPP) == Message{PointMass}
    @test isApplicable(SPMultiplicationOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPMultiplicationOutNPP(nothing, Message(Univariate, PointMass, m=2.0), Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPMultiplicationOutNPP(nothing, Message(Multivariate, PointMass, m=[2.0]), Message(MatrixVariate, PointMass, m=mat(1.0))) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPMultiplicationIn1GNP" begin
    @test SPMultiplicationIn1GNP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1GNP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPMultiplicationIn1GNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPMultiplicationIn1GNP(Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=3.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=12.0)
    @test ruleSPMultiplicationIn1GNP(Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(3.0)), nothing, Message(Multivariate, PointMass, m=[2.0])) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=12.0 + tiny)
    @test ruleSPMultiplicationIn1GNP(Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(3.0)), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[2.0], w=mat(12.0 + tiny))
end

@testset "SPMultiplicationIn1PNP" begin
    @test SPMultiplicationIn1PNP <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationIn1PNP) == Message{PointMass}
    @test isApplicable(SPMultiplicationIn1PNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPMultiplicationIn1PNP(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, PointMass, m=2.0)) == Message(Univariate, PointMass, m=0.5)
    @test ruleSPMultiplicationIn1PNP(Message(Multivariate, PointMass, m=[1.0]), nothing, Message(MatrixVariate, PointMass, m=mat(2.0))) == Message(Multivariate, PointMass, m=[0.5])
end

@testset "SPMultiplicationAGPN" begin
    @test SPMultiplicationAGPN <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationAGPN) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPMultiplicationAGPN, [Message{Gaussian}, Message{PointMass}, Nothing])

    @test ruleSPMultiplicationAGPN(Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=3.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=12.0)
end

@testset "SPMultiplicationAPPN" begin
    @test SPMultiplicationAPPN <: SumProductRule{Multiplication}
    @test outboundType(SPMultiplicationAPPN) == Message{PointMass}
    @test isApplicable(SPMultiplicationAPPN, [Message{PointMass}, Message{PointMass}, Nothing])

    @test ruleSPMultiplicationAPPN(Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=2.0), nothing) == Message(Univariate, PointMass, m=0.5)
end

end # module
