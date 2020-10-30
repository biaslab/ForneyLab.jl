module ProbitTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, requiresBreaker, breakerParameters
using ForneyLab: SPProbitOutNG, SPProbitIn1PN, EPProbitIn1BG, EPProbitIn1CG, EPProbitIn1PG
using StatsFuns: normcdf

@testset "requiresBreaker and breakerParameters" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    ng = GaussianMeanVariance(x, 0.0, 1.0)
    np = Probit(y, x)
    
    @test !requiresBreaker(np.i[:out]) # Dangling
    @test !requiresBreaker(np.i[:in1])
    @test !requiresBreaker(ng.i[:m])
    @test !requiresBreaker(ng.i[:out]) # requiresBreaker not applicable for EP

    @test_throws Exception breakerParameters(np.i[:out])
    @test breakerParameters(ng.i[:out]) == (Message{GaussianMeanVariance, Univariate}, ())
end


#-------------
# Update rules
#-------------

@testset "SPProbitOutNG" begin
    @test SPProbitOutNG <: SumProductRule{Probit}
    @test outboundType(SPProbitOutNG) == Message{Bernoulli}
    @test isApplicable(SPProbitOutNG, [Nothing, Message{Gaussian}]) 
    @test !isApplicable(SPProbitOutNG, [Message{Bernoulli}, Nothing])

    @test ruleSPProbitOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, Bernoulli, p=normcdf(1/sqrt(1+0.5)))
end

@testset "SPProbitIn1PN" begin
    @test SPProbitIn1PN <: SumProductRule{Probit}
    @test outboundType(SPProbitIn1PN) == Message{Function}
    @test isApplicable(SPProbitIn1PN, [Message{PointMass}, Nothing]) 

    @test isa(ruleSPProbitIn1PN(Message(Univariate, PointMass, m=1.0), nothing), Message{Function, Univariate})
end

@testset "EPProbitIn1BG" begin
    @test EPProbitIn1BG <: ExpectationPropagationRule{Probit}
    @test outboundType(EPProbitIn1BG) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPProbitIn1BG, [Message{Bernoulli}, Message{Gaussian}], 2) 
    @test !isApplicable(EPProbitIn1BG, [Message{PointMass}, Message{Gaussian}], 2)

    @test ruleEPProbitIn1BG(Message(Univariate, Bernoulli, p=1.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPProbitIn1BG(Message(Univariate, Bernoulli, p=0.8), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPProbitIn1BG(Message(Univariate, Bernoulli, p=0.5), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "EPProbitIn1CG" begin
    @test EPProbitIn1CG <: ExpectationPropagationRule{Probit}
    @test outboundType(EPProbitIn1CG) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPProbitIn1CG, [Message{Categorical}, Message{Gaussian}], 2) 
    @test !isApplicable(EPProbitIn1CG, [Message{PointMass}, Message{Gaussian}], 2)

    @test ruleEPProbitIn1CG(Message(Univariate, Categorical, p=[1.0, 0.0]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPProbitIn1CG(Message(Univariate, Categorical, p=[0.8, 0.2]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPProbitIn1CG(Message(Univariate, Categorical, p=[0.5, 0.5]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "EPProbitIn1PG" begin
    @test EPProbitIn1PG <: ExpectationPropagationRule{Probit}
    @test outboundType(EPProbitIn1PG) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPProbitIn1PG, [Message{PointMass}, Message{Gaussian}], 2) 
    @test !isApplicable(EPProbitIn1PG, [Message{Bernoulli}, Message{Gaussian}], 2) 

    @test ruleEPProbitIn1PG(Message(Univariate, PointMass, m=true), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPProbitIn1PG(Message(Univariate, PointMass, m=NaN), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "averageEnergy" begin
    @test averageEnergy(Probit, ProbabilityDistribution(Univariate, Bernoulli, p=1.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == -1.291942516803335
    @test averageEnergy(Probit, ProbabilityDistribution(Univariate, PointMass, m=1.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == -1.291942516803335
end

end # module