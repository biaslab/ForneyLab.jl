module ProbitTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, requiresBreaker, breakerParameters
using ForneyLab: SPProbitOutNG, EPProbitIn1GB, EPProbitIn1GC, EPProbitIn1GP

@testset "requiresBreaker and breakerParameters" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    ng = GaussianMeanVariance(x, 0.0, 1.0)
    np = Probit(y, x)
    
    @test !requiresBreaker(np.i[:out]) # Dangling
    @test !requiresBreaker(np.i[:in1])
    @test !requiresBreaker(ng.i[:m])
    @test requiresBreaker(ng.i[:out])

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

    @test ruleSPProbitOutNG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, Bernoulli, p=ForneyLab.Φ(1/sqrt(1+0.5)))
end

@testset "EPProbitIn1GB" begin
    @test EPProbitIn1GB <: ExpectationPropagationRule{Probit}
    @test outboundType(EPProbitIn1GB) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPProbitIn1GB, [Message{Bernoulli}, Message{Gaussian}], 2) 
    @test !isApplicable(EPProbitIn1GB, [Message{PointMass}, Message{Gaussian}], 2)

    @test ruleEPProbitIn1GB(Message(Univariate, Bernoulli, p=1.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPProbitIn1GB(Message(Univariate, Bernoulli, p=0.8), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPProbitIn1GB(Message(Univariate, Bernoulli, p=0.5), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "EPProbitIn1GC" begin
    @test EPProbitIn1GC <: ExpectationPropagationRule{Probit}
    @test outboundType(EPProbitIn1GC) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPProbitIn1GC, [Message{Categorical}, Message{Gaussian}], 2) 
    @test !isApplicable(EPProbitIn1GC, [Message{PointMass}, Message{Gaussian}], 2)

    @test ruleEPProbitIn1GC(Message(Univariate, Categorical, p=[1.0, 0.0]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPProbitIn1GC(Message(Univariate, Categorical, p=[0.8, 0.2]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPProbitIn1GC(Message(Univariate, Categorical, p=[0.5, 0.5]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "EPProbitIn1GP" begin
    @test EPProbitIn1GP <: ExpectationPropagationRule{Probit}
    @test outboundType(EPProbitIn1GP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPProbitIn1GP, [Message{PointMass}, Message{Gaussian}], 2) 
    @test !isApplicable(EPProbitIn1GP, [Message{Bernoulli}, Message{Gaussian}], 2) 

    @test ruleEPProbitIn1GP(Message(Univariate, PointMass, m=true), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPProbitIn1GP(Message(Univariate, PointMass, m=NaN), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "averageEnergy" begin
    @test averageEnergy(Probit, ProbabilityDistribution(Univariate, Bernoulli, p=1.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == -1.291942516803335
    @test averageEnergy(Probit, ProbabilityDistribution(Univariate, PointMass, m=1.0), ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == -1.291942516803335
end

end # module