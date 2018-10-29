module SigmoidTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, mapToBernoulliParameterRange
import ForneyLab: SPSigmoidBinVG, EPSigmoidRealGB, EPSigmoidRealGC, EPSigmoidRealGP

@testset "mapToBernoulliParameterRange" begin
    @test mapToBernoulliParameterRange(1.0) == 1.0
    @test mapToBernoulliParameterRange(0.5) == 0.75
    @test mapToBernoulliParameterRange(0.0) == 0.5
    @test mapToBernoulliParameterRange(-1.0) == 0.0
    @test mapToBernoulliParameterRange(true) == 1.0
    @test mapToBernoulliParameterRange(false) == 0.0
    @test mapToBernoulliParameterRange(NaN) == 0.5
    @test_throws Exception mapToBernoulliParameterRange(2.0)
end


#-------------
# Update rules
#-------------

@testset "SPSigmoidBinVG" begin
    @test SPSigmoidBinVG <: SumProductRule{Sigmoid}
    @test outboundType(SPSigmoidBinVG) == Message{Bernoulli}
    @test isApplicable(SPSigmoidBinVG, [Nothing, Message{Gaussian}]) 
    @test !isApplicable(SPSigmoidBinVG, [Message{Bernoulli}, Nothing])

    @test ruleSPSigmoidBinVG(nothing, Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, Bernoulli, p=ForneyLab.Φ(1/sqrt(1+0.5)))
end

@testset "EPSigmoidRealGB" begin
    @test EPSigmoidRealGB <: ExpectationPropagationRule{Sigmoid}
    @test outboundType(EPSigmoidRealGB) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPSigmoidRealGB, [Message{Bernoulli}, Message{Gaussian}], 2) 
    @test !isApplicable(EPSigmoidRealGB, [Message{PointMass}, Message{Gaussian}], 2)

    @test ruleEPSigmoidRealGB(Message(Univariate, Bernoulli, p=1.0), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPSigmoidRealGB(Message(Univariate, Bernoulli, p=0.8), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPSigmoidRealGB(Message(Univariate, Bernoulli, p=0.5), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "EPSigmoidRealGC" begin
    @test EPSigmoidRealGC <: ExpectationPropagationRule{Sigmoid}
    @test outboundType(EPSigmoidRealGC) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPSigmoidRealGC, [Message{Categorical}, Message{Gaussian}], 2) 
    @test !isApplicable(EPSigmoidRealGC, [Message{PointMass}, Message{Gaussian}], 2)

    @test ruleEPSigmoidRealGC(Message(Univariate, Categorical, p=[1.0, 0.0]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPSigmoidRealGC(Message(Univariate, Categorical, p=[0.8, 0.2]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPSigmoidRealGC(Message(Univariate, Categorical, p=[0.5, 0.5]), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

@testset "EPSigmoidRealGP" begin
    @test EPSigmoidRealGP <: ExpectationPropagationRule{Sigmoid}
    @test outboundType(EPSigmoidRealGP) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(EPSigmoidRealGP, [Message{PointMass}, Message{Gaussian}], 2) 
    @test !isApplicable(EPSigmoidRealGP, [Message{Bernoulli}, Message{Gaussian}], 2) 

    @test ruleEPSigmoidRealGP(Message(Univariate, PointMass, m=true), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPSigmoidRealGP(Message(Univariate, PointMass, m=NaN), Message(Univariate, GaussianMeanVariance, m=1.0, v=0.5)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4e-12, w=4e-12)
end

end # module