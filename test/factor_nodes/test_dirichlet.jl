module DirichletTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeVar, vague, dims
using ForneyLab: SPDirichletOutNP, VBDirichletOut, VBDirichletIn1
using SpecialFunctions: digamma

@testset "Dirichlet ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0]) == ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[2.0, 3.0]))
    @test ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0]) == ProbabilityDistribution{MatrixVariate, Dirichlet}(Dict(:a=>[2.0 3.0; 4.0 5.0]))
    @test_throws Exception ProbabilityDistribution(Univariate, Dirichlet)
    @test ProbabilityDistribution(Dirichlet, a=[2.0, 3.0]) == ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[2.0, 3.0]))
    @test ProbabilityDistribution(Dirichlet) == ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0]))
    @test Message(Dirichlet) == Message{Dirichlet, Multivariate}(ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0])))
    @test Message(Multivariate, Dirichlet) == Message{Dirichlet, Multivariate}(ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0])))
    @test Message(MatrixVariate, Dirichlet) == Message{Dirichlet, MatrixVariate}(ProbabilityDistribution{MatrixVariate, Dirichlet}(Dict(:a=>ones(3,3))))
    @test_throws Exception Message(Univariate, Dirichlet)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0, 2.0])) == 3
    @test dims(ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 2.0 2.0; 2.0 2.0 2.0])) == (2,3)
end

@testset "vague" begin
    @test vague(Dirichlet, 3) == ProbabilityDistribution(Dirichlet, a=ones(3))
    @test vague(Dirichlet, (3,)) == ProbabilityDistribution(Dirichlet, a=ones(3))
    @test vague(Dirichlet, (2,3)) == ProbabilityDistribution(MatrixVariate, Dirichlet, a=ones(2,3))
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0])) == [0.5, 0.5]
    @test unsafeMean(ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 1.0; 2.0 3.0])) == [0.5 0.25; 0.5 0.75]
    @test unsafeLogMean(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0])) == [digamma(2.0), digamma(3.0)] .- digamma(5.0)
    @test unsafeLogMean(ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 4.0; 3.0 5.0])) == [digamma(2.0) digamma(4.0); digamma(3.0) digamma(5.0)] - [digamma(5.0) digamma(9.0); digamma(5.0) digamma(9.0)]
    @test unsafeVar(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0])) == [0.05, 0.05]
end

@testset "prod!" begin
    # Multivariate
    @test ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0]) == ProbabilityDistribution(Multivariate, Dirichlet, a=[3.0, 4.0])
    @test ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0]) * ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]) == ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])
    @test ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0]) == ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])
    @test_throws Exception ProbabilityDistribution(Multivariate, PointMass, m=[-0.1, 0.8, 0.1]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0])
    @test_throws Exception ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.9, 0.1]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0])

    # MatrixVariate
    @test ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 2.0; 3.0 3.0]) * ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0]) == ProbabilityDistribution(MatrixVariate, Dirichlet, a=[3.0 4.0; 6.0 7.0])
    @test ProbabilityDistribution(MatrixVariate, Dirichlet, a=[1.0 3.0; 2.0 4.0]) * ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7]) == ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7])
    @test ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7]) * ProbabilityDistribution(MatrixVariate, Dirichlet, a=[1.0 3.0; 2.0 4.0]) == ProbabilityDistribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7])
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(Multivariate, Dirichlet, a=[0.2,3.0,1.5]),[2,3,7]), 3.2556382883760024)
    @test isapprox(logPdf(ProbabilityDistribution(MatrixVariate, Dirichlet, a=[0.2 1.4; 3.0 1.8]),[2 7; 3 3]), 3.0442561618507087)
end

#-------------
# Update rules
#-------------

@testset "SPDirichletOutNP" begin
    @test SPDirichletOutNP <: SumProductRule{Dirichlet}
    @test outboundType(SPDirichletOutNP) == Message{Dirichlet}
    @test isApplicable(SPDirichletOutNP, [Nothing, Message{PointMass}])

    @test ruleSPDirichletOutNP(nothing, Message(Multivariate, PointMass, m=[2.0, 3.0])) == Message(Multivariate, Dirichlet, a=[2.0, 3.0])
    @test ruleSPDirichletOutNP(nothing, Message(MatrixVariate, PointMass, m=[2.0 3.0; 4.0 5.0])) == Message(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0])
end

@testset "VBDirichletOut" begin
    @test VBDirichletOut <: NaiveVariationalRule{Dirichlet}
    @test outboundType(VBDirichletOut) == Message{Dirichlet}
    @test isApplicable(VBDirichletOut, [Nothing, ProbabilityDistribution])

    @test ruleVBDirichletOut(nothing, ProbabilityDistribution(Multivariate, PointMass, m=[2.0, 3.0])) == Message(Multivariate, Dirichlet, a=[2.0, 3.0])
    @test ruleVBDirichletOut(nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=[2.0 3.0; 4.0 5.0])) == Message(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0])
end

@testset "VBDirichletIn1" begin
    @test VBDirichletIn1 <: NaiveVariationalRule{Dirichlet}
    @test outboundType(VBDirichletIn1) == Message{Function}
    @test isApplicable(VBDirichletIn1, [ProbabilityDistribution, Nothing])
end

@testset "averageEnergy and differentialEntropy" begin
    # Multivariate
    @test isapprox(differentialEntropy(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0])), averageEnergy(Dirichlet, ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0]), ProbabilityDistribution(Multivariate, PointMass, m=[2.0, 3.0])))
    @test isapprox(averageEnergy(Dirichlet, ProbabilityDistribution(Multivariate, Dirichlet, a=[4.0, 5.0]), ProbabilityDistribution(Multivariate, PointMass, m=[2.0, 3.0])), averageEnergy(Beta, ProbabilityDistribution(Univariate, Beta, a=4.0, b=5.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)))
    @test isapprox(averageEnergy(Dirichlet,ProbabilityDistribution(Multivariate,Dirichlet,a=[2.0,3.0,4.0]),ProbabilityDistribution(Multivariate,SampleList,s=[[1.0,2.0,1.0],[3.,3.,1.],[2.,2.,2]],w=[0.1,0.7,0.3])),0.5093876870003795)

    # MatrixVariate
    @test differentialEntropy(ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0])) == differentialEntropy(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 4.0])) + differentialEntropy(ProbabilityDistribution(Multivariate, Dirichlet, a=[3.0, 5.0]))
    @test averageEnergy(Dirichlet, ProbabilityDistribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0]), ProbabilityDistribution(MatrixVariate, PointMass, m=[6.0 7.0; 8.0 9.0])) == averageEnergy(Dirichlet, ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 4.0]), ProbabilityDistribution(Multivariate, PointMass, m=[6.0, 8.0])) + averageEnergy(Dirichlet, ProbabilityDistribution(Multivariate, Dirichlet, a=[3.0, 5.0]), ProbabilityDistribution(Multivariate, PointMass, m=[7.0, 9.0]))
end

end # module
