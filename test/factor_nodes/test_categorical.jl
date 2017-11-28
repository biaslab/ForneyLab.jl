module CategoricalTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, vague, dims
import ForneyLab: SPCategoricalOutVP, VBCategoricalOut

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

@testset "prod!" begin
    @test ProbabilityDistribution(Categorical, p=[0.2, 0.8])*ProbabilityDistribution(Categorical, p=[0.8, 0.2]) == ProbabilityDistribution(Categorical, p=[0.5, 0.5])    
    @test ProbabilityDistribution(Categorical, p=[0.25, 0.5, 0.25]) * ProbabilityDistribution(Categorical, p=[1/3, 1/3, 1/3]) == ProbabilityDistribution(Categorical, p=[0.25, 0.5, 0.25])
    @test_throws Exception ProbabilityDistribution(Categorical, p=[0.0, 0.5, 0.5]) * ProbabilityDistribution(Categorical, p=[1.0, 0.0, 0.0])
end

#-------------
# Update rules
#-------------

@testset "SPCategoricalOutVP" begin
    @test SPCategoricalOutVP <: SumProductRule{Categorical}
    @test outboundType(SPCategoricalOutVP) == Message{Categorical}
    @test isApplicable(SPCategoricalOutVP, [Void, Message{PointMass}]) 

    @test ruleSPCategoricalOutVP(nothing, Message(Multivariate, PointMass, m=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[0.1, 0.8, 0.1])
end

@testset "VBCategoricalOut" begin
    @test VBCategoricalOut <: VariationalRule{Categorical}
    @test outboundType(VBCategoricalOut) == Message{Categorical}
    @test isApplicable(VBCategoricalOut, [Void, ProbabilityDistribution])
    @test !isApplicable(VBCategoricalOut, [ProbabilityDistribution, Void])

    @test ruleVBCategoricalOut(nothing, ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[0.10000000000000002, 0.8, 0.10000000000000002])
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.8, 0.1])) == averageEnergy(Categorical, ProbabilityDistribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]), ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]))
    @test differentialEntropy(ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8])) == differentialEntropy(ProbabilityDistribution(Univariate, Bernoulli, p=0.2))
    @test averageEnergy(Categorical, ProbabilityDistribution(Univariate, Categorical, p=[0.2, 0.8]), ProbabilityDistribution(Multivariate, PointMass, m=[0.1, 0.9])) == averageEnergy(Bernoulli, ProbabilityDistribution(Univariate, Bernoulli, p=0.2), ProbabilityDistribution(Univariate, PointMass, m=0.1))
end

end # module