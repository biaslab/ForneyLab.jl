module SigmoidTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, mapToBernoulliParameterRange, SPSigmoidGV, EPSigmoidGB1, EPSigmoidGP1

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

@testset "SPSigmoidGV" begin
    @test SPSigmoidGV <: SumProductRule{Sigmoid}
    @test outboundType(SPSigmoidGV) == Message{Bernoulli}
    @test isApplicable(SPSigmoidGV, [Message{Gaussian}, Void]) 
    @test !isApplicable(SPSigmoidGV, [Void, Message{Bernoulli}])

    @test ruleSPSigmoidGV(Message(Gaussian, m=1.0, v=0.5), nothing) == Message(Bernoulli, p=ForneyLab.Φ(1/sqrt(1+0.5)))
end

@testset "EPSigmoidGB1" begin
    @test EPSigmoidGB1 <: ExpectationPropagationRule{Sigmoid}
    @test outboundType(EPSigmoidGB1) == Message{Gaussian}
    @test isApplicable(EPSigmoidGB1, [Message{Gaussian}, Message{Bernoulli}], 1) 
    @test !isApplicable(EPSigmoidGB1, [Message{Gaussian}, Message{PointMass}], 1) 

    @test ruleEPSigmoidGB1(Message(Gaussian, m=1.0, v=0.5), Message(Bernoulli, p=1.0)) == Message(Gaussian, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPSigmoidGB1(Message(Gaussian, m=1.0, v=0.5), Message(Bernoulli, p=0.8)) == Message(Gaussian, xi=0.4270174959448596, w=0.19914199922339604)
    @test ruleEPSigmoidGB1(Message(Gaussian, m=1.0, v=0.5), Message(Bernoulli, p=0.5)) == Message(Gaussian, xi=0.0, w=1e-12)
end

@testset "EPSigmoidGP1" begin
    @test EPSigmoidGP1 <: ExpectationPropagationRule{Sigmoid}
    @test outboundType(EPSigmoidGP1) == Message{Gaussian}
    @test isApplicable(EPSigmoidGP1, [Message{Gaussian}, Message{PointMass}], 1) 
    @test !isApplicable(EPSigmoidGP1, [Message{Gaussian}, Message{Bernoulli}], 1) 

    @test ruleEPSigmoidGP1(Message(Gaussian, m=1.0, v=0.5), Message(PointMass, m=true)) == Message(Gaussian, xi=0.6723616582693994, w=0.3295003993960708)
    @test ruleEPSigmoidGP1(Message(Gaussian, m=1.0, v=0.5), Message(PointMass, m=NaN)) == Message(Gaussian, xi=0.0, w=1e-12)
end

end # module