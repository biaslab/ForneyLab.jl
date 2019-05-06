module AdditionTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPAdditionOutNGG, SPAdditionOutNGP, SPAdditionOutNPG, SPAdditionIn1GNG, SPAdditionIn1PNG, SPAdditionIn2GGN, SPAdditionIn2PGN, SPAdditionIn1GNP, SPAdditionIn2GPN, SPAdditionOutNPP, SPAdditionIn1PNP, SPAdditionIn2PPN

@testset "Addition node construction through + syntax" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV y ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x + y
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)

    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x + 1.0
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
end

@testset "Addition node construction through - syntax" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV y ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x - y
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
    # This syntax works by changing the order interfaces are attached in. Below checks whether the pairings are succesful by
    # matching out and in2 interfaces of the corresponding nodes
    @test g.nodes[:addition_1].i[:out] == g.nodes[:gaussianmeanvariance_1].i[:out].partner
    @test g.nodes[:addition_1].i[:in2] == g.nodes[:gaussianmeanvariance_2].i[:out].partner

    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV z = x - 1.0
    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
end

#-------------
# Update rules
#-------------

@testset "SPAdditionOutNGG" begin
    @test SPAdditionOutNGG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutNGG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionOutNGG, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test !isApplicable(SPAdditionOutNGG, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPAdditionOutNGG, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleSPAdditionOutNGG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0)) == Message(Univariate, GaussianMeanVariance, m=4.0, v=6.0)
    @test ruleSPAdditionOutNGG(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0))) == Message(Multivariate, GaussianMeanVariance, m=[4.0], v=mat(6.0))
end

@testset "SPAdditionIn2GGN" begin
    @test SPAdditionIn2GGN <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GGN) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn2GGN, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleSPAdditionIn2GGN(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.0, v=6.0)
    @test ruleSPAdditionIn2GGN(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(6.0))
end

@testset "SPAdditionIn1GNG" begin
    @test SPAdditionIn1GNG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn1GNG, [Message{Gaussian}, Nothing, Message{Gaussian}])

    @test ruleSPAdditionIn1GNG(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=6.0)
    @test ruleSPAdditionIn1GNG(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(6.0))
end

@testset "SPAdditionOutNGP" begin
    @test SPAdditionOutNGP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutNGP) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionOutNGP, [Nothing, Message{Gaussian}, Message{PointMass}])

    @test ruleSPAdditionOutNGP(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
    @test ruleSPAdditionOutNGP(nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), Message(Multivariate, PointMass, m=[3.0])) == Message(Multivariate, GaussianMeanVariance, m=[4.0], v=mat(2.0))
end

@testset "SPAdditionOutNPG" begin
    @test SPAdditionOutNPG <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutNPG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionOutNPG, [Nothing, Message{PointMass}, Message{Gaussian}])

    @test ruleSPAdditionOutNPG(nothing, Message(Univariate, PointMass, m=3.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, GaussianMeanVariance, m=4.0, v=2.0)
    @test ruleSPAdditionOutNPG(nothing, Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[4.0], v=mat(2.0))
end

@testset "SPAdditionIn1PNG" begin
    @test SPAdditionIn1PNG <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PNG) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn1PNG, [Message{PointMass}, Nothing, Message{Gaussian}])

    @test ruleSPAdditionIn1PNG(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0)
    @test ruleSPAdditionIn1PNG(Message(Multivariate, PointMass, m=[3.0]), nothing, Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0))) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(2.0))
end

@testset "SPAdditionIn2PGN" begin
    @test SPAdditionIn2PGN <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PGN) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn2PGN, [Message{PointMass}, Message{Gaussian}, Nothing])

    @test ruleSPAdditionIn2PGN(Message(Univariate, PointMass, m=3.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=2.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.0, v=2.0)
    @test ruleSPAdditionIn2PGN(Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, GaussianMeanVariance, m=[1.0], v=mat(2.0)), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(2.0))
end

@testset "SPAdditionIn1GNP" begin
    @test SPAdditionIn1GNP <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1GNP) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn1GNP, [Message{Gaussian}, Nothing, Message{PointMass}])

    @test ruleSPAdditionIn1GNP(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, GaussianMeanVariance, m=2.0, v=4.0)
    @test ruleSPAdditionIn1GNP(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), nothing, Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))
end

@testset "SPAdditionIn2GPN" begin
    @test SPAdditionIn2GPN <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2GPN) == Message{GaussianMeanVariance}
    @test isApplicable(SPAdditionIn2GPN, [Message{Gaussian}, Message{PointMass}, Nothing])

    @test ruleSPAdditionIn2GPN(Message(Univariate, GaussianMeanVariance, m=3.0, v=4.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, GaussianMeanVariance, m=2.0, v=4.0)
    @test ruleSPAdditionIn2GPN(Message(Multivariate, GaussianMeanVariance, m=[3.0], v=mat(4.0)), Message(Multivariate, PointMass, m=[1.0]), nothing) == Message(Multivariate, GaussianMeanVariance, m=[2.0], v=mat(4.0))
end

@testset "SPAdditionOutNPP" begin
    @test SPAdditionOutNPP <: SumProductRule{Addition}
    @test outboundType(SPAdditionOutNPP) == Message{PointMass}
    @test isApplicable(SPAdditionOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPAdditionOutNPP(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, PointMass, m=3.0)) == Message(Univariate, PointMass, m=4.0)
    @test ruleSPAdditionOutNPP(nothing, Message(Multivariate, PointMass, m=[1.0]), Message(Multivariate, PointMass, m=[3.0])) == Message(Multivariate, PointMass, m=[4.0])
end

@testset "SPAdditionIn2PPN" begin
    @test SPAdditionIn2PPN <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn2PPN) == Message{PointMass}
    @test isApplicable(SPAdditionIn2PPN, [Message{PointMass}, Message{PointMass}, Nothing])

    @test ruleSPAdditionIn2PPN(Message(Univariate, PointMass, m=3.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPAdditionIn2PPN(Message(Multivariate, PointMass, m=[3.0]), Message(Multivariate, PointMass, m=[1.0]), nothing) == Message(Multivariate, PointMass, m=[2.0])
end

@testset "SPAdditionIn1PNP" begin
    @test SPAdditionIn1PNP <: SumProductRule{Addition}
    @test outboundType(SPAdditionIn1PNP) == Message{PointMass}
    @test isApplicable(SPAdditionIn1PNP, [Message{PointMass}, Nothing, Message{PointMass}])

    @test ruleSPAdditionIn1PNP(Message(Univariate, PointMass, m=3.0), nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, PointMass, m=2.0)
    @test ruleSPAdditionIn1PNP(Message(Multivariate, PointMass, m=[3.0]), nothing, Message(Multivariate, PointMass, m=[1.0])) == Message(Multivariate, PointMass, m=[2.0])
end

end # module
