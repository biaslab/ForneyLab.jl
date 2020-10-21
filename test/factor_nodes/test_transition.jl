module TransitionTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPTransitionOutNPP, SPTransitionIn1PNP, SPTransitionOutNCP, SPTransitionIn1CNP, VBTransitionOut, VBTransitionIn1, VBTransitionA, SVBTransitionOutVCD, SVBTransitionIn1CVD, SVBTransitionADV, MTransitionCCD, MTransitionCCN


#-------------
# Update rules
#-------------

@testset "SPTransitionOutNPP" begin
    @test SPTransitionOutNPP <: SumProductRule{Transition}
    @test outboundType(SPTransitionOutNPP) == Message{Categorical}
    @test isApplicable(SPTransitionOutNPP, [Nothing, Message{PointMass}, Message{PointMass}]) 
    @test !isApplicable(SPTransitionOutNPP, [Message{PointMass}, Nothing, Message{PointMass}]) 

    @test ruleSPTransitionOutNPP(nothing, Message(Multivariate, PointMass, m=[0.1, 0.4, 0.5]), Message(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.3660714285714285, 0.27678571428571425, 0.35714285714285715])
end

@testset "SPTransitionIn1PNP" begin
    @test SPTransitionIn1PNP <: SumProductRule{Transition}
    @test outboundType(SPTransitionIn1PNP) == Message{Categorical}
    @test isApplicable(SPTransitionIn1PNP, [Message{PointMass}, Nothing, Message{PointMass}]) 

    @test ruleSPTransitionIn1PNP(Message(Multivariate, PointMass, m=[0.1, 0.4, 0.5]), nothing, Message(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.23000000000000004, 0.43, 0.33999999999999997])
end

@testset "SPTransitionOutNCP" begin
    @test SPTransitionOutNCP <: SumProductRule{Transition}
    @test outboundType(SPTransitionOutNCP) == Message{Categorical}
    @test isApplicable(SPTransitionOutNCP, [Nothing, Message{Categorical}, Message{PointMass}]) 
    @test !isApplicable(SPTransitionOutNCP, [Message{Categorical}, Nothing, Message{PointMass}]) 

    @test ruleSPTransitionOutNCP(nothing, Message(Univariate, Categorical, p=[0.1, 0.4, 0.5]), Message(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.3660714285714285, 0.27678571428571425, 0.35714285714285715])
end

@testset "SPTransitionIn1CNP" begin
    @test SPTransitionIn1CNP <: SumProductRule{Transition}
    @test outboundType(SPTransitionIn1CNP) == Message{Categorical}
    @test isApplicable(SPTransitionIn1CNP, [Message{Categorical}, Nothing, Message{PointMass}]) 

    @test ruleSPTransitionIn1CNP(Message(Univariate, Categorical, p=[0.1, 0.4, 0.5]), nothing, Message(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.23000000000000004, 0.43, 0.33999999999999997])
end

@testset "VBTransitionOut" begin
    @test VBTransitionOut <: NaiveVariationalRule{Transition}
    @test outboundType(VBTransitionOut) == Message{Categorical}
    @test isApplicable(VBTransitionOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution]) 

    @test ruleVBTransitionOut(nothing, ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.4, 0.5]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.2994385327821292, 0.3260399327754116, 0.37452153444245917])
end

@testset "VBTransitionIn1" begin
    @test VBTransitionIn1 <: NaiveVariationalRule{Transition}
    @test outboundType(VBTransitionIn1) == Message{Categorical}
    @test isApplicable(VBTransitionIn1, [ProbabilityDistribution, Nothing, ProbabilityDistribution]) 

    @test ruleVBTransitionIn1(ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.4, 0.5]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.208905960977562, 0.4255474215521219, 0.365546617470316])
end

@testset "VBTransitionA" begin
    @test VBTransitionA <: NaiveVariationalRule{Transition}
    @test outboundType(VBTransitionA) == Message{Dirichlet}
    @test isApplicable(VBTransitionA, [ProbabilityDistribution, ProbabilityDistribution, Nothing]) 

    @test ruleVBTransitionA(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.7, 0.1]), ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.4, 0.5]), nothing) == Message(MatrixVariate, Dirichlet, a=[1.02 1.08 1.1; 1.07 1.28 1.35; 1.01 1.04 1.05])
end

@testset "SVBTransitionOutVCD" begin
    @test SVBTransitionOutVCD <: StructuredVariationalRule{Transition}
    @test outboundType(SVBTransitionOutVCD) == Message{Categorical}
    @test isApplicable(SVBTransitionOutVCD, [Nothing, Message{Categorical}, ProbabilityDistribution]) 

    @test ruleSVBTransitionOutVCD(nothing, Message(Univariate, Categorical, p=[0.1, 0.4, 0.5]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.3660714285714286, 0.27678571428571425, 0.35714285714285715])
end

@testset "SVBTransitionIn1CVD" begin
    @test SVBTransitionIn1CVD <: StructuredVariationalRule{Transition}
    @test outboundType(SVBTransitionIn1CVD) == Message{Categorical}
    @test isApplicable(SVBTransitionIn1CVD, [Message{Categorical}, Nothing, ProbabilityDistribution]) 

    @test ruleSVBTransitionIn1CVD(Message(Univariate, Categorical, p=[0.1, 0.4, 0.5]), nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == Message(Univariate, Categorical, p=[0.23000000000000007, 0.42999999999999994, 0.33999999999999997])
end

@testset "SVBTransitionADV" begin
    @test SVBTransitionADV <: StructuredVariationalRule{Transition}
    @test outboundType(SVBTransitionADV) == Message{Dirichlet}
    @test isApplicable(SVBTransitionADV, [ProbabilityDistribution, Nothing]) 

    @test ruleSVBTransitionADV(ProbabilityDistribution(Multivariate, Contingency, p=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3]), nothing) == Message(MatrixVariate, Dirichlet, a=[1.2 1.1 1.7; 1.4 1.3 1.3; 1.1 1.6 1.3])
end

@testset "MTransitionCCD" begin
    @test MTransitionCCD <: MarginalRule{Transition}
    @test isApplicable(MTransitionCCD, [Message{Categorical}, Message{Categorical}, ProbabilityDistribution]) 
    @test !isApplicable(MTransitionCCD, [Message{Categorical}, Message{Categorical}, Nothing]) 

    @test ruleMTransitionCCD(Message(Univariate, Categorical, p=[0.2, 0.7, 0.1]), Message(Univariate, Categorical, p=[0.1, 0.4, 0.5]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == ProbabilityDistribution(Multivariate, Contingency, p=[0.011799410029498528 0.023598820058997057 0.20648967551622419; 0.08259587020648967 0.247787610619469 0.3097345132743362; 0.002949852507374632 0.07079646017699116 0.04424778761061947])
end

@testset "MTransitionCCN" begin
    @test MTransitionCCN <: MarginalRule{Transition}
    @test isApplicable(MTransitionCCN, [Message{Categorical}, Message{Categorical}, Nothing]) 
    @test !isApplicable(MTransitionCCN, [Message{Categorical}, Message{Categorical}, ProbabilityDistribution]) 

    @test ruleMTransitionCCN(Message(Univariate, Categorical, p=[0.2, 0.7, 0.1]), Message(Univariate, Categorical, p=[0.1, 0.4, 0.5]), Message(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == ProbabilityDistribution(Multivariate, Contingency, p=[0.011799410029498528 0.023598820058997057 0.20648967551622419; 0.08259587020648967 0.247787610619469 0.3097345132743362; 0.002949852507374632 0.07079646017699116 0.04424778761061947])
end

@testset "averageEnergy" begin
    @test averageEnergy(Transition, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.7, 0.1]), ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.4, 0.5]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == 1.1783637941354863
    @test averageEnergy(Transition, ProbabilityDistribution(Multivariate, Contingency, p=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3]), ProbabilityDistribution(MatrixVariate, PointMass, m=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == 2.7886642527453405
end

end