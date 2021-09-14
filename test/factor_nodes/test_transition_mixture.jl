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

    @test ruleSPTransitionMixtureOutNCPX(nothing, Message(Univariate, Categorical, p=[0.1, 0.9]), Message(Multivariate, PointMass, m=[1.0, 0.0]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.19000000000000003, 0.81])
end

@testset "SPTransitionMixtureIn1CNPX" begin
    @test SPTransitionMixtureIn1CNPX <: SumProductRule{TransitionMixture}
    @test outboundType(SPTransitionMixtureIn1CNPX) == Message{Categorical}
    @test isApplicable(SPTransitionMixtureIn1CNPX, [Message{Categorical}, Nothing, Message{PointMass}, Message{PointMass}, Message{PointMass}]) 

    @test ruleSPTransitionMixtureIn1CNPX(Message(Univariate, Categorical, p=[0.2, 0.8]), nothing, Message(Multivariate, PointMass, m=[1.0, 0.0]), Message(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), Message(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5211267605633803, 0.4788732394366197])
end

@testset "VBTransitionMixtureZ" begin
    @test VBTransitionMixtureZ <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureZ) == Message{Categorical}
    @test isApplicable(VBTransitionMixtureZ, [ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureZ(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.9]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5201135171653145, 0.47988648283468543])
end

@testset "VBTransitionMixtureOut" begin
    @test VBTransitionMixtureOut <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureOut) == Message{Categorical}
    @test isApplicable(VBTransitionMixtureOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureOut(nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.9]), ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.319739310361466, 0.6802606896385339])
end

@testset "VBTransitionMixtureIn1" begin
    @test VBTransitionMixtureIn1 <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureIn1) == Message{Categorical}
    @test isApplicable(VBTransitionMixtureIn1, [ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureIn1(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(Univariate, Categorical, p=[0.5081812668233174, 0.4918187331766825])
end

@testset "VBTransitionMixtureA" begin
    @test VBTransitionMixtureA <: NaiveVariationalRule{TransitionMixture}
    @test outboundType(VBTransitionMixtureA) == Message{Dirichlet}
    @test isApplicable(VBTransitionMixtureA, [ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution]) 

    @test ruleVBTransitionMixtureA(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.9]), ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == Message(MatrixVariate, Dirichlet, a=[1.006 1.054; 1.024 1.216])
end

@testset "averageEnergy" begin
    @test averageEnergy(TransitionMixture, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.9]), ProbabilityDistribution(Univariate, Categorical, p=[0.3, 0.7]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.2; 0.9 0.8]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.3 0.4; 0.7 0.6])) == 0.5611909800043979
end

end #module