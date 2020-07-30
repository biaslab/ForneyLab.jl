module CategoricalTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, vague, dims, logPdf
using ForneyLab: SPCategoricalOutNP, VBCategoricalOut, VBCategoricalIn1
using SparseArrays: SparseVector, spzeros

@testset "Categorical ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]) == ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>[0.1, 0.8, 0.1]))
    @test_throws Exception ProbabilityDistribution(Multivariate, Categorical)
    @test ProbabilityDistribution(Categorical, p=[0.1, 0.8, 0.1]) == ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>[0.1, 0.8, 0.1]))
    @test ProbabilityDistribution(Categorical) == ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>[1/3, 1/3, 1/3]))
    @test Message(Categorical) == Message{Categorical, Univariate}(ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>[1/3, 1/3, 1/3])))
    @test Message(Univariate, Categorical) == Message{Categorical, Univariate}(ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>[1/3, 1/3, 1/3])))
    @test_throws Exception Message(Multivariate, Categorical)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Categorical, p=[0.1, 0.8, 0.1])) == 1
end

@testset "vague" begin
    @test vague(Categorical) == ProbabilityDistribution(Categorical, p=[1/3, 1/3, 1/3])
    @test vague(Categorical, 4) == ProbabilityDistribution(Categorical, p=[1/4, 1/4, 1/4, 1/4])
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Categorical, p=[0.1, 0.8, 0.1])) == [0.1, 0.8, 0.1]
end

@testset "sample" begin
    s = sample(ProbabilityDistribution(Categorical, p=[0.1, 0.8, 0.1]))
    @test isa(s, SparseVector)
    @test sum(s.==1.0) == 1
    @test sum(s.==0.0) == 2
end

@testset "log pdf" begin
    x = spzeros(Float64, 3)
    x[3] = 1.0
    @test isapprox(logPdf(ProbabilityDistribution(Categorical, p=[0.2, 0.3, 0.5]),x), -0.6931471805599453)
end

@testset "prod!" begin
    @test ProbabilityDistribution(Categorical, p=[0.2, 0.8])*ProbabilityDistribution(Categorical, p=[0.8, 0.2]) == ProbabilityDistribution(Categorical, p=[0.5, 0.5])
    @test ProbabilityDistribution(Categorical, p=[0.25, 0.5, 0.25]) * ProbabilityDistribution(Categorical, p=[1/3, 1/3, 1/3]) == ProbabilityDistribution(Categorical, p=[0.25, 0.5, 0.25])
    @test ProbabilityDistribution(Categorical, p=[0.0, 0.5, 0.5]) * ProbabilityDistribution(Categorical, p=[1.0, 0.0, 0.0]) == ProbabilityDistribution(Categorical, p=ones(3)/3)
end

#-------------
# Update rules
#-------------

@testset "SPCategoricalOutNP" begin
    @test SPCategoricalOutNP <: SumProductRule{Categorical}
    @test outboundType(SPCategoricalOutNP) == Message{Categorical}
    @test isApplicable(SPCategoricalOutNP, [Nothing, Message{PointMass}])

    @test ruleSPCategoricalOutNP(nothing, Message(Multivariate, PointMass, m=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[0.1, 0.8, 0.1])
end

@testset "VBCategoricalOut" begin
    @test VBCategoricalOut <: NaiveVariationalRule{Categorical}
    @test outboundType(VBCategoricalOut) == Message{Categorical}
    @test isApplicable(VBCategoricalOut, [Nothing, ProbabilityDistribution])
    @test !isApplicable(VBCategoricalOut, [ProbabilityDistribution, Nothing])

    @test ruleVBCategoricalOut(nothing, ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[0.10000000000000002, 0.8, 0.10000000000000002])
    @test ruleVBCategoricalOut(nothing, ProbabilityDistribution(Multivariate, Dirichlet, a=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[7.799215056092699e-5, 0.999844015698878, 7.799215056092699e-5])
end

@testset "VBCategoricalIn1" begin
    @test VBCategoricalIn1 <: NaiveVariationalRule{Categorical}
    @test outboundType(VBCategoricalIn1) == Message{Dirichlet}
    @test isApplicable(VBCategoricalIn1, [ProbabilityDistribution, Nothing])

    @test ruleVBCategoricalIn1(ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]), nothing) == Message(Multivariate, Dirichlet, a=[1.1, 1.8, 1.1])
    @test ruleVBCategoricalIn1(ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]), nothing) == Message(Multivariate, Dirichlet, a=[1.1, 1.8, 1.1])
    @test ruleVBCategoricalIn1(ProbabilityDistribution(Univariate, Bernoulli, p=0.1), nothing) == Message(Multivariate, Dirichlet, a=[1.1, 1.9])
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.8, 0.1])) == averageEnergy(Categorical, ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]), ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]))
    @test differentialEntropy(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8])) == differentialEntropy(ProbabilityDistribution(Univariate, Bernoulli, p=0.2))
    @test averageEnergy(Categorical, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.9])) == averageEnergy(Bernoulli, ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, PointMass, m=0.1))
end

end # module
