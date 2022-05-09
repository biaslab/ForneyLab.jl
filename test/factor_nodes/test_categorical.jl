module CategoricalTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, vague, dims, logPdf, naturalParams, standardDistribution
using ForneyLab: SPCategoricalOutNP, VBCategoricalOut, VBCategoricalIn1
using SparseArrays: SparseVector, spzeros

@testset "Categorical Distribution and Message construction" begin
    @test Distribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]) == Distribution{Univariate, Categorical}(Dict(:p=>[0.1, 0.8, 0.1]))
    @test_throws Exception Distribution(Multivariate, Categorical)
    @test Distribution(Categorical, p=[0.1, 0.8, 0.1]) == Distribution{Univariate, Categorical}(Dict(:p=>[0.1, 0.8, 0.1]))
    @test Distribution(Categorical) == Distribution{Univariate, Categorical}(Dict(:p=>[1/3, 1/3, 1/3]))
    @test Message(Categorical) == Message{Categorical, Univariate}(Distribution{Univariate, Categorical}(Dict(:p=>[1/3, 1/3, 1/3])))
    @test Message(Univariate, Categorical) == Message{Categorical, Univariate}(Distribution{Univariate, Categorical}(Dict(:p=>[1/3, 1/3, 1/3])))
    @test_throws Exception Message(Multivariate, Categorical)
end

@testset "dims" begin
    @test dims(Distribution(Categorical, p=[0.1, 0.8, 0.1])) == ()
end

@testset "vague" begin
    @test vague(Categorical) == Distribution(Categorical, p=[1/3, 1/3, 1/3])
    @test vague(Categorical, 4) == Distribution(Categorical, p=[1/4, 1/4, 1/4, 1/4])
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Distribution(Categorical, p=[0.1, 0.8, 0.1])) == [0.1, 0.8, 0.1]
end

@testset "sample" begin
    s = sample(Distribution(Categorical, p=[0.1, 0.8, 0.1]))
    @test isa(s, SparseVector)
    @test sum(s.==1.0) == 1
    @test sum(s.==0.0) == 2
end

@testset "log pdf" begin
    x = spzeros(Float64, 3)
    x[3] = 1.0
    @test isapprox(logPdf(Distribution(Categorical, p=[0.2, 0.3, 0.5]),x), -0.6931471805599453)
end

@testset "prod!" begin
    @test Distribution(Categorical, p=[0.2, 0.8])*Distribution(Categorical, p=[0.8, 0.2]) == Distribution(Categorical, p=[0.5, 0.5])
    @test Distribution(Categorical, p=[0.2, 0.8])*Distribution(Bernoulli, p=0.8) == Distribution(Categorical, p=[0.5, 0.5])
    @test Distribution(Categorical, p=[0.25, 0.5, 0.25]) * Distribution(Categorical, p=[1/3, 1/3, 1/3]) == Distribution(Categorical, p=[0.25, 0.5, 0.25])
    @test Distribution(Categorical, p=[0.0, 0.5, 0.5]) * Distribution(Categorical, p=[1.0, 0.0, 0.0]) == Distribution(Categorical, p=ones(3)/3)
end

@testset "natural parameters" begin
    d = Distribution(Univariate, Categorical, p=[0.2, 0.8])
    η = naturalParams(d)
    s = standardDistribution(Univariate, Categorical, η=η)
    @test d.params[:p] == s.params[:p] # Test conversion consistency

    x = [[0.0, 1.0], [1.0, 0.0]]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, Categorical, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
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
    @test isApplicable(VBCategoricalOut, [Nothing, Distribution])
    @test !isApplicable(VBCategoricalOut, [Distribution, Nothing])

    @test ruleVBCategoricalOut(nothing, Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[0.10000000000000002, 0.8, 0.10000000000000002])
    @test ruleVBCategoricalOut(nothing, Distribution(Multivariate, Dirichlet, a=[0.1, 0.8, 0.1])) == Message(Univariate, Categorical, p=[7.799215056092699e-5, 0.999844015698878, 7.799215056092699e-5])
end

@testset "VBCategoricalIn1" begin
    @test VBCategoricalIn1 <: NaiveVariationalRule{Categorical}
    @test outboundType(VBCategoricalIn1) == Message{Dirichlet}
    @test isApplicable(VBCategoricalIn1, [Distribution, Nothing])

    @test ruleVBCategoricalIn1(Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]), nothing) == Message(Multivariate, Dirichlet, a=[1.1, 1.8, 1.1])
    @test ruleVBCategoricalIn1(Distribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]), nothing) == Message(Multivariate, Dirichlet, a=[1.1, 1.8, 1.1])
    @test ruleVBCategoricalIn1(Distribution(Univariate, Bernoulli, p=0.1), nothing) == Message(Multivariate, Dirichlet, a=[1.1, 1.9])
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, Categorical, p=[0.1, 0.8, 0.1])) == averageEnergy(Categorical, Distribution(Univariate, Categorical, p=[0.1, 0.8, 0.1]), Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]))
    @test differentialEntropy(Distribution(Univariate, Categorical, p=[0.2, 0.8])) == differentialEntropy(Distribution(Univariate, Bernoulli, p=0.2))
    @test averageEnergy(Categorical, Distribution(Univariate, Categorical, p=[0.2, 0.8]), Distribution(Multivariate, PointMass, m=[0.1, 0.9])) == averageEnergy(Bernoulli, Distribution(Univariate, Bernoulli, p=0.2), Distribution(Univariate, PointMass, m=0.1))
end

end # module
