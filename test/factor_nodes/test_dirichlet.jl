module DirichletTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeVar, vague, dims
import ForneyLab: SPDirichletOutVP, VBDirichletOut

@testset "Dirichlet ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0]) == ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[2.0, 3.0]))
    @test_throws Exception ProbabilityDistribution(Univariate, Dirichlet)
    @test ProbabilityDistribution(Dirichlet, a=[2.0, 3.0]) == ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[2.0, 3.0]))
    @test ProbabilityDistribution(Dirichlet) == ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0]))
    @test Message(Dirichlet) == Message{Dirichlet, Multivariate}(ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0])))
    @test Message(Multivariate, Dirichlet) == Message{Dirichlet, Multivariate}(ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0])))
    @test_throws Exception Message(Univariate, Dirichlet)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Dirichlet, a=[2.0, 2.0, 2.0])) == 3
end

@testset "vague" begin
    @test vague(Dirichlet, 3) == ProbabilityDistribution(Dirichlet, a=ones(3))
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0])) == [0.5, 0.5]
    @test unsafeLogMean(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0])) == [digamma(2.0), digamma(3.0)] - digamma(5.0)
    @test unsafeVar(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0])) == [0.05, 0.05]
end

@testset "prod!" begin
    @test ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 2.0]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0]) == ProbabilityDistribution(Multivariate, Dirichlet, a=[3.0, 4.0])
    @test ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0]) * ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]) == ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])
    @test ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0]) == ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])
    @test_throws Exception ProbabilityDistribution(Multivariate, PointMass, m=[-0.1, 0.8, 0.1]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0])
    @test_throws Exception ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.9, 0.1]) * ProbabilityDistribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0])
end

#-------------
# Update rules
#-------------

@testset "SPDirichletOutVP" begin
    @test SPDirichletOutVP <: SumProductRule{Dirichlet}
    @test outboundType(SPDirichletOutVP) == Message{Dirichlet}
    @test isApplicable(SPDirichletOutVP, [Void, Message{PointMass}])

    @test ruleSPDirichletOutVP(nothing, Message(Multivariate, PointMass, m=[2.0, 3.0])) == Message(Multivariate, Dirichlet, a=[2.0, 3.0])
end

@testset "VBDirichletOut" begin
    @test VBDirichletOut <: NaiveVariationalRule{Dirichlet}
    @test outboundType(VBDirichletOut) == Message{Dirichlet}
    @test isApplicable(VBDirichletOut, [Void, ProbabilityDistribution])

    @test ruleVBDirichletOut(nothing, ProbabilityDistribution(Multivariate, PointMass, m=[2.0, 3.0])) == Message(Multivariate, Dirichlet, a=[2.0, 3.0])
end

@testset "averageEnergy and differentialEntropy" begin
    @test isapprox(differentialEntropy(ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0])), averageEnergy(Dirichlet, ProbabilityDistribution(Multivariate, Dirichlet, a=[2.0, 3.0]), ProbabilityDistribution(Multivariate, PointMass, m=[2.0, 3.0])))
    @test isapprox(averageEnergy(Dirichlet, ProbabilityDistribution(Multivariate, Dirichlet, a=[4.0, 5.0]), ProbabilityDistribution(Multivariate, PointMass, m=[2.0, 3.0])), averageEnergy(Beta, ProbabilityDistribution(Univariate, Beta, a=4.0, b=5.0), ProbabilityDistribution(Univariate, PointMass, m=2.0), ProbabilityDistribution(Univariate, PointMass, m=3.0)))
end

end # module