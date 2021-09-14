module TransitionMixtureTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPTransitionMixtureOutNCPX, SPTransitionMixtureIn1CNPX, VBTransitionMixtureZ, VBTransitionMixtureOut, VBTransitionMixtureIn1, VBTransitionMixtureA

#-------------
# Update rules
#-------------

@testset "SPTransitionMixtureOutNCPX" begin
    @test SPTransitionMixtureOutNCPX <: SumProductRule{TransitionMixture}
    @test outboundType(SPTransitionMixtureOutNCPX) == Message{Categorical}
    @test isApplicable(SPTransitionMixtureOutNCPX, [Nothing, Message{Categorical}, Message{PointMass}, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPTransitionMixtureOutNCPX(nothing, Message(Univariate, Categorical, p=[0.1, 0.9]), Message(Univariate, PointMass, p=[1.0, 0.0]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.11, 0.89])
end

@testset "SPTransitionMixtureIn1NCPX" begin
    @test SPTransitionMixtureIn1NCPX <: SumProductRule{TransitionMixture}
    @test outboundType(SPTransitionMixtureIn1NCPX) == Message{Categorical}
    @test isApplicable(SPTransitionMixtureIn1NCPX, [Message{Categorical}, Nothing, Message{PointMass}, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPTransitionMixtureIn1NCPX(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, Message(Univariate, PointMass, p=[1.0, 0.0]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.11, 0.89])
end

@testset "VBTransitionMixtureZ" begin
    @test VBTransitionMixtureZ <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureZ) == Message{Categorical}
    @test isApplicable(VBTransitionMixtureZ, [ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureZ(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), Message(Univariate, Categorical, p=[0.1, 0.9]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.11, 0.89])
end

@testset "VBTransitionMixtureOut" begin
    @test VBTransitionMixtureOut <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureOut) == Message{Categorical}
    @test isApplicable(VBTransitionMixtureOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureOut(nothing, Message(Univariate, Categorical, p=[0.1, 0.9]), ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.11, 0.89])
end

@testset "VBTransitionMixtureIn1" begin
    @test VBTransitionMixtureIn1 <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureIn1) == Message{Categorical}
    @test isApplicable(VBTransitionMixtureIn1, [ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureIn1(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.11, 0.89])
end

@testset "VBTransitionMixtureA" begin
    @test VBTransitionMixtureA <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureA) == Message{Dirichlet}
    @test isApplicable(VBTransitionMixtureA, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureA(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), Message(Univariate, Categorical, p=[0.1, 0.9]), ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(MatrixVariate, Dirichlet, a=[1.0 2.0; 9.0 8.0])
end

@testset "averageEnergy" begin
    @test averageEnergy(TransitionMixture, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.9]), ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == 1.0
end