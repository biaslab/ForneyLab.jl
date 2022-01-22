module BernoulliTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar, vague, dims, logPdf, naturalParams, standardDistribution
using ForneyLab: SPBernoulliOutNP, SPBernoulliIn1PN, SPBernoulliOutNB, VBBernoulliOut, VBBernoulliIn1

@testset "Bernoulli Distribution and Message construction" begin
    @test Distribution(Univariate, Bernoulli, p=0.2) == Distribution{Univariate, Bernoulli}(Dict(:p=>0.2))
    @test_throws Exception Distribution(Multivariate, Bernoulli)
    @test Distribution(Bernoulli, p=0.2) == Distribution{Univariate, Bernoulli}(Dict(:p=>0.2))
    @test Distribution(Bernoulli) == Distribution{Univariate, Bernoulli}(Dict(:p=>0.5))
    @test Message(Bernoulli) == Message{Bernoulli, Univariate}(Distribution{Univariate, Bernoulli}(Dict(:p=>0.5)))
    @test Message(Univariate, Bernoulli) == Message{Bernoulli, Univariate}(Distribution{Univariate, Bernoulli}(Dict(:p=>0.5)))
    @test_throws Exception Message(Multivariate, Bernoulli)
end

@testset "dims" begin
    @test dims(Distribution(Bernoulli, p=0.5)) == ()
end

@testset "vague" begin
    @test vague(Bernoulli) == Distribution(Bernoulli, p=0.5)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Distribution(Bernoulli, p=0.2)) == 0.2
    @test unsafeVar(Distribution(Bernoulli, p=0.5)) == 0.25
end

@testset "log pmf" begin
    @test isapprox(logPdf(Distribution(Bernoulli, p=0.2),1), -1.6094379124341003)
end

@testset "prod!" begin
    @test Distribution(Bernoulli, p=0.2) * Distribution(Bernoulli, p=0.8) == Distribution(Bernoulli, p=0.5000000000000001)
    @test_throws Exception Distribution(Bernoulli, p=0.0) * Distribution(Bernoulli, p=1.0)
end

@testset "natural parameters" begin
    d = Distribution(Univariate, Bernoulli, p=0.2)
    η = naturalParams(d)
    s = standardDistribution(Univariate, Bernoulli, η=η)
    @test d.params[:p] == s.params[:p] # Test conversion consistency
    
    x = [0.0, 1.0]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, Bernoulli, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
end


#-------------
# Update rules
#-------------

@testset "SPBernoulliOutNP" begin
    @test SPBernoulliOutNP <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliOutNP) == Message{Bernoulli}
    @test isApplicable(SPBernoulliOutNP, [Nothing, Message{PointMass}])

    @test ruleSPBernoulliOutNP(nothing, Message(Univariate, PointMass, m=1.0)) == Message(Univariate, Bernoulli, p=1.0)
end

@testset "SPBernoulliIn1PN" begin
    @test SPBernoulliIn1PN <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliIn1PN) == Message{Beta}
    @test isApplicable(SPBernoulliIn1PN, [Message{PointMass}, Nothing])

    @test ruleSPBernoulliIn1PN(Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, Beta, a=2.0, b=1.0)
end

@testset "SPBernoulliOutNB" begin
    @test SPBernoulliOutNB <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliOutNB) == Message{Bernoulli}
    @test isApplicable(SPBernoulliOutNB, [Nothing, Message{Beta}])

    @test ruleSPBernoulliOutNB(nothing, Message(Univariate, Beta, a=1.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.5)
end

@testset "VBBernoulliOut" begin
    @test VBBernoulliOut <: NaiveVariationalRule{Bernoulli}


    @test outboundType(VBBernoulliOut) == Message{Bernoulli}
    @test isApplicable(VBBernoulliOut, [Nothing, Distribution])
    @test !isApplicable(VBBernoulliOut, [Distribution, Nothing])

    @test ruleVBBernoulliOut(nothing, Distribution(Univariate, PointMass, m=0.2)) == Message(Univariate, Bernoulli, p=0.2)
    @test ruleVBBernoulliOut(nothing, Distribution(Univariate, Beta, a=1.0, b=1.0)) == Message(Univariate, Bernoulli, p=0.5)
end

@testset "VBBernoulliIn1" begin
    @test VBBernoulliIn1 <: NaiveVariationalRule{Bernoulli}
    @test outboundType(VBBernoulliIn1) == Message{Beta}
    @test isApplicable(VBBernoulliIn1, [Distribution, Nothing])

    @test ruleVBBernoulliIn1(Distribution(Univariate, Bernoulli, p=0.2), nothing) == Message(Univariate, Beta, a=1.2, b=1.8)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(Univariate, Bernoulli, p=0.25)) == averageEnergy(Bernoulli, Distribution(Univariate, Bernoulli, p=0.25), Distribution(Univariate, PointMass, m=0.25))
end

end # module
