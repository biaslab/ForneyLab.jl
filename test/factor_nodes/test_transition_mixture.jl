module TransitionMixtureTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPTransitionMixtureOutNCCPX, SPTransitionMixtureIn1CNCPX, SPTransitionMixtureZCCNPX, SVBTransitionMixtureOutNCCDX, SVBTransitionMixtureIn1CNCDX, SVBTransitionMixtureZCCNDX,  SVBTransitionMixtureA, MTransitionMixtureCCCDX, MTransitionMixtureCCCNX

#-------------
# Update rules
#-------------

@testset "SPTransitionMixtureOutNCCPX" begin
    @test SPTransitionMixtureOutNCCPX <: SumProductRule{TransitionMixture}
    @test outboundType(SPTransitionMixtureOutNCCPX) == Message{Categorical}
    @test isApplicable(SPTransitionMixtureOutNCCPX, [Nothing, Message{Categorical}, Message{Categorical}, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPTransitionMixtureOutNCCPX(nothing, Message(Univariate, Categorical, p=[0.1, 0.9]), Message(Univariate, Categorical, p=[1.0, 0.0]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.19000000000000003, 0.81])
end

@testset "SPTransitionMixtureIn1CNCPX" begin
    @test SPTransitionMixtureIn1CNCPX <: SumProductRule{TransitionMixture}
    @test outboundType(SPTransitionMixtureIn1CNCPX) == Message{Categorical}
    @test isApplicable(SPTransitionMixtureIn1CNCPX, [Message{Categorical}, Nothing, Message{Categorical}, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPTransitionMixtureIn1CNCPX(Message(Univariate, Categorical, p=[0.2, 0.8]), nothing, Message(Univariate, Categorical, p=[1.0, 0.0]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5211267605633803, 0.4788732394366197])
end

@testset "SPTransitionMixtureZCCNPX" begin
    @test SPTransitionMixtureZCCNPX <: SumProductRule{TransitionMixture}
    @test outboundType(SPTransitionMixtureZCCNPX) == Message{Categorical}
    @test isApplicable(SPTransitionMixtureZCCNPX, [Message{Categorical}, Message{Categorical}, Nothing, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPTransitionMixtureZCCNPX(Message(Univariate, Categorical, p=[0.2, 0.8]), Message(Univariate, Categorical, p=[1.0, 0.0]), nothing, Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5441176470588237, 0.4558823529411764])
end

@testset "SVBTransitionMixtureOutNCCDX" begin
    @test SVBTransitionMixtureOutNCCDX <: StructuredVariationalRule{TransitionMixture}
    @test outboundType(SVBTransitionMixtureOutNCCDX) == Message{Categorical}
    @test isApplicable(SVBTransitionMixtureOutNCCDX, [Nothing, Message{Categorical}, Message{Categorical}, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleSVBTransitionMixtureOutNCCDX(nothing, Message(Univariate, Categorical, p=[0.1, 0.9]), Message(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.33000000000033997, 0.6699999999996599])
end

@testset "SVBTransitionMixtureIn1CNCDX" begin
    @test SVBTransitionMixtureIn1CNCDX <: StructuredVariationalRule{TransitionMixture}
    @test outboundType(SVBTransitionMixtureIn1CNCDX) == Message{Categorical}
    @test isApplicable(SVBTransitionMixtureIn1CNCDX, [Message{Categorical}, Nothing, Message{Categorical}, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleSVBTransitionMixtureIn1CNCDX(Message(Univariate, Categorical, p=[0.2, 0.8]), nothing, Message(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5239616613418148, 0.4760383386581853])
end

@testset "SVBTransitionMixtureZCCNDX" begin
    @test SVBTransitionMixtureZCCNDX <: StructuredVariationalRule{TransitionMixture}
    @test outboundType(SVBTransitionMixtureZCCNDX) == Message{Categorical}
    @test isApplicable(SVBTransitionMixtureZCCNDX, [Message{Categorical}, Message{Categorical}, Nothing, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleSVBTransitionMixtureZCCNDX(Message(Univariate, Categorical, p=[0.2, 0.8]), Message(Univariate, Categorical, p=[0.1, 0.9]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5479233226837061, 0.4520766773162939])
end

@testset "SVBTransitionMixtureA" begin
    @test SVBTransitionMixtureA <: StructuredVariationalRule{TransitionMixture}
    @test outboundType(SVBTransitionMixtureA) == Message{Dirichlet}
    @test isApplicable(SVBTransitionMixtureA, [ProbabilityDistribution, Nothing, ProbabilityDistribution]) 

    @test ruleSVBTransitionMixtureA(ProbabilityDistribution(Multivariate, Contingency, p=[0.3*[0.1 0.6; 0.2 0.1], 0.7*[0.1 0.6; 0.2 0.1]]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(MatrixVariate, Dirichlet, a=[1.03 1.18; 1.06 1.03])
end

@testset "MTransitionMixtureCCCDX" begin
    @test MTransitionMixtureCCCDX <: MarginalRule{TransitionMixture}
    @test isApplicable(MTransitionMixtureCCCDX, [Message{Categorical}, Message{Categorical}, Message{Categorical}, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleMTransitionMixtureCCCDX(Message(Univariate, Categorical, p=[0.2, 0.8]), Message(Univariate, Categorical, p=[0.1, 0.9]), Message(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == ProbabilityDistribution(Multivariate, Contingency, p=[[0.0009966777408637877 0.017940199335548173; 0.035880398671096346 0.28704318936877077], [0.00697674418604651 0.08372093023255814; 0.06511627906976744 0.5023255813953489]])
end

@testset "MTransitionMixtureCCCNX" begin
    @test MTransitionMixtureCCCNX <: MarginalRule{TransitionMixture}
    @test isApplicable(MTransitionMixtureCCCNX, [Message{Categorical}, Message{Categorical}, Message{Categorical}, Nothing, Nothing]) 

    @test ruleMTransitionMixtureCCCNX(Message(Univariate, Categorical, p=[0.2, 0.8]), Message(Univariate, Categorical, p=[0.1, 0.9]), Message(Univariate, Categorical, p=[0.3, 0.7]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == ProbabilityDistribution(Multivariate, Contingency, p=[[0.0009966777408637875 0.017940199335548173; 0.035880398671096346 0.28704318936877077], [0.0069767441860465115 0.08372093023255814; 0.06511627906976744 0.5023255813953489]])
end

@testset "averageEnergy" begin
    @test averageEnergy(TransitionMixture, ProbabilityDistribution(Multivariate, Contingency, p=[0.3*[0.1 0.6; 0.2 0.1], 0.7*[0.1 0.6; 0.2 0.1]]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == 0.9266048040118577
end

end #module