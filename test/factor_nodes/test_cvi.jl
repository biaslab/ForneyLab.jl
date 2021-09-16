module CviTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar
using ForneyLab: SPCVIOutNFactorNode, SPCVIIn1Factor, SPCVIOutNFactorNodeX, SPCVIInFactorX, MCVIFactorX
using Flux.Optimise

f1(x) = x
f2(x::Real,z::Real) = x*z
f2(x::Vector,z::Vector) = x.*z
f3(x) = [x,2*x]
f4(x,z,w) = x*z*w
opt = Descent(0.1)

isUnivariateSetSampleList(dist::ProbabilityDistribution{Univariate, SetSampleList}) = true
isUnivariateSetSampleList(dist::ProbabilityDistribution{Multivariate, SetSampleList}) = false
isUnivariateSampleList(dist::ProbabilityDistribution{Univariate, SampleList}) = true
isUnivariateGaussian(dist::ProbabilityDistribution{Univariate, F}) where F<:Gaussian = true
isMultivariateGaussian(dist::ProbabilityDistribution{Multivariate, F}) where F<:Gaussian = true
isGamma(dist::ProbabilityDistribution{Univariate, Gamma}) = true
isBernoulli(dist::ProbabilityDistribution{Univariate, Bernoulli}) = true

@testset "Cvi Node Construction" begin
    fg = FactorGraph()
    x = Variable()
    z = Variable()
    y = Variable()
    @test Cvi(y, x, g=f1, opt=opt, num_samples=10, num_iterations=10).q_memory == [ProbabilityDistribution(Univariate,GaussianMeanVariance,m=0,v=1)]
    @test Cvi(y, x, g=f1, opt=opt, num_samples=10, num_iterations=10).online_inference == false
    @test Cvi(y, x, z, g=f2, opt=opt, num_samples=10, num_iterations=[10,100], q=[vague(GaussianMeanVariance),vague(GaussianMeanPrecision)]).q_memory == [vague(GaussianMeanVariance),vague(GaussianMeanPrecision)]
    @test Cvi(y, x, z, g=f2, opt=opt, num_samples=10, num_iterations=[10,100], q=[vague(GaussianMeanVariance),vague(GaussianMeanPrecision)]).online_inference == [false,false]
    @test Cvi(y, x, z, g=f2, opt=opt, num_samples=10, num_iterations=[10,1], q=[vague(GaussianMeanVariance),vague(GaussianMeanPrecision)],online_inference=[false,true]).online_inference == [false,true]
end

#-------------
# Update rules
#-------------

@testset "SPCVIOutNFactorNode" begin
    @test SPCVIOutNFactorNode <: SumProductRule{CVI}
    @test outboundType(SPCVIOutNFactorNode) == Message{SetSampleList}
    @test isApplicable(SPCVIOutNFactorNode, [Nothing, Message{Gaussian}])
    @test isApplicable(SPCVIOutNFactorNode, [Nothing, Message{Bernoulli}])

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    node1 =  Cvi(y, x, g=f1, opt=opt, num_samples=1000, num_iterations=10)
    node2 =  Cvi(y, x, g=f3, opt=opt, num_samples=1000, num_iterations=10)
    msg1 = ruleSPCVIOutNFactorNode(node1.id,nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny))
    msg2 = ruleSPCVIOutNFactorNode(node2.id,nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny))

    @test isUnivariateSetSampleList(msg1.dist)
    @test !isUnivariateSetSampleList(msg2.dist)

    @test msg1.dist.params[:node_id] == node1.id
    @test msg2.dist.params[:node_id] == node2.id
end

@testset "SPCVIIn1Factor" begin
    @test SPCVIIn1Factor <: SumProductRule{CVI}
    @test outboundType(SPCVIIn1Factor) == Message{FactorNode}
    @test isApplicable(SPCVIIn1Factor, [Message{Gaussian}, Nothing])
    @test isApplicable(SPCVIIn1Factor, [Message{Bernoulli}, Nothing])
    @test isApplicable(SPCVIIn1Factor, [Message{FactorFunction}, Nothing])

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    node1 =  Cvi(y, x, g=f1, opt=opt, num_samples=1000, num_iterations=1)
    node2 =  Cvi(y, x, g=f3, opt=opt, num_samples=1000, num_iterations=1)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=0., v=10.)
    msg_in2 = Message(Univariate, Gamma, a=1., b=1.)
    msg_in3 = Message(Multivariate, GaussianMeanPrecision, m=zeros(2), w=[1. 0.;0. 0.1])
    msg1 = ruleSPCVIIn1Factor(node1.id, Message(Univariate, GaussianMeanPrecision, m=2.0, w=1.5), msg_in1)
    msg2 = ruleSPCVIIn1Factor(node2.id, Message(Multivariate, GaussianMeanVariance, m=zeros(2), v=[1. 0.;0. 1.]), msg_in1)
    msg3 = ruleSPCVIIn1Factor(node2.id, Message(Multivariate, GaussianMeanVariance, m=zeros(2), v=[1. 0.;0. 1.]), msg_in2)
    msg4 = ruleSPCVIIn1Factor(node1.id, Message(Multivariate, GaussianMeanVariance, m=[2.0,-1.], v=[1. 0.;0. 1.]), msg_in3)

    @test isUnivariateGaussian(msg1.dist)
    @test isUnivariateGaussian(msg2.dist)
    @test isGamma(msg3.dist)
    @test isMultivariateGaussian(msg4.dist)
    @test isMultivariateGaussian(node1.q[1])
    @test isMultivariateGaussian(node1.q_memory[1])
end

@testset "SPCVIOutNFactorNodeX" begin
    @test SPCVIOutNFactorNodeX <: SumProductRule{CVI}
    @test outboundType(SPCVIOutNFactorNodeX) == Message{SetSampleList}
    @test isApplicable(SPCVIOutNFactorNodeX, [Nothing, Message{Gaussian}, Message{Gaussian}])
    @test isApplicable(SPCVIOutNFactorNodeX, [Nothing, Message{Bernoulli}, Message{Gamma}])

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    y = Variable()
    node1 =  Cvi(y, x, z, g=f2, opt=opt, num_samples=1000, num_iterations=10)
    msg1 = ruleSPCVIOutNFactorNodeX(node1.id, nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny), Message(Univariate, Gamma, a=1., b=1.))
    msg2 = ruleSPCVIOutNFactorNodeX(node1.id, nothing, Message(Multivariate, GaussianMeanPrecision, m=[-1.2,2.0], w=[0.5 0.2;0.2 0.4]), Message(Multivariate, GaussianMeanPrecision, m=[-1.2,2.0], w=[0.5 0.2;0.2 0.4]))

    @test isUnivariateSetSampleList(msg1.dist)
    @test !isUnivariateSetSampleList(msg2.dist)

    @test msg1.dist.params[:node_id] == node1.id
    @test msg2.dist.params[:node_id] == node1.id
end

@testset "SPCVIInFactorX" begin
    @test SPCVIInFactorX <: SumProductRule{CVI}
    @test outboundType(SPCVIInFactorX) == Message{FactorNode}
    @test isApplicable(SPCVIInFactorX, [Message{Gaussian}, Message{Gaussian}, Nothing])
    @test isApplicable(SPCVIInFactorX, [Message{Gaussian}, Message{Bernoulli}, Message{Gamma}, Nothing])

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    y = Variable()
    node1 =  Cvi(y, x, z, g=f2, opt=[opt,opt], num_samples=1000, num_iterations=[1,1], q=[vague(GaussianMeanVariance), vague(Gamma)])
    msg_out = Message(Univariate, GaussianMeanPrecision, m=1.5, w=0.2)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny)
    msg_in2 = Message(Univariate, Gamma, a=1., b=1.)

    msg1 = ruleSPCVIInFactorX(node1.id, 1, msg_out, msg_in1, msg_in2)
    @test isUnivariateGaussian(msg1.dist)
    @test node1.q[1] != vague(GaussianMeanVariance)
    @test node1.q[2] != vague(Gamma)
    new_q = [node1.q[1], node1.q[2]]
    @test node1.q_memory[1] != vague(GaussianMeanVariance)
    @test node1.q_memory[2] != vague(Gamma)
    @test node1.infer_memory == 1
    msg2 = ruleSPCVIInFactorX(node1.id, 2, msg_out, msg_in1, msg_in2)
    @test isGamma(msg2.dist)
    @test node1.infer_memory == 0
    @test node1.q[1] == new_q[1]
    @test node1.q[2] == new_q[2]

    node2 =  Cvi(y, x, z, g=f2, opt=[opt,opt], num_samples=1000, num_iterations=[1,1], q=[vague(GaussianMeanVariance), vague(Gamma)], online_inference=[true,false])
    msg1 = ruleSPCVIInFactorX(node2.id, 1, msg_out, msg_in1, msg_in2)
    @test isUnivariateGaussian(msg1.dist)
    @test node2.q[1] != vague(GaussianMeanVariance)
    @test node2.q[2] != vague(Gamma)
    @test node2.q_memory[1] == vague(GaussianMeanVariance)
    @test node2.q_memory[2] != vague(Gamma)
    @test node2.infer_memory == 1
    msg2 = ruleSPCVIInFactorX(node2.id, 2, msg_out, msg_in1, msg_in2)
    @test isGamma(msg2.dist)
    @test node2.infer_memory == 0

    w = Variable()
    node3 =  Cvi(y, x, z, w, g=f4, opt=[opt,opt,opt], num_samples=1000, num_iterations=[1,1,1], q=[vague(GaussianMeanVariance), vague(Gamma), vague(Bernoulli)], online_inference=[false,true,false])
    msg_in3 = Message(Univariate, Bernoulli, p=.9)
    msg1 = ruleSPCVIInFactorX(node3.id, 1, msg_out, msg_in1, msg_in2, msg_in3)
    @test isUnivariateGaussian(msg1.dist)
    @test node3.q[1] != vague(GaussianMeanVariance)
    @test node3.q[2] != vague(Gamma)
    @test node3.q[3] != vague(Bernoulli)
    @test node3.q_memory[1] != vague(GaussianMeanVariance)
    @test node3.q_memory[2] == vague(Gamma)
    @test node3.q_memory[3] != vague(Bernoulli)
    @test node3.infer_memory == 2
    msg2 = ruleSPCVIInFactorX(node3.id, 2, msg_out, msg_in1, msg_in2, msg_in3)
    @test isGamma(msg2.dist)
    @test node3.infer_memory == 1
    msg3 = ruleSPCVIInFactorX(node3.id, 3, msg_out, msg_in1, msg_in2, msg_in3)
    @test isBernoulli(msg3.dist)
    @test node3.infer_memory == 0
end

#-------------
# SetSampleList
#-------------

@testset "SetSampleList" begin
    fg = FactorGraph()
    x = Variable()
    z = Variable()
    w = Variable()
    y = Variable()
    node =  Cvi(y, x, z, w, g=f4, opt=[opt,opt,opt], num_samples=1000, num_iterations=[1,1,1], q=[vague(GaussianMeanVariance), vague(Gamma), vague(Bernoulli)], online_inference=[false,true,false])
    msg_out = Message(Univariate, GaussianMeanPrecision, m=1.5, w=0.2)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny)
    msg_in2 = Message(Univariate, Gamma, a=1., b=1.)
    msg_in3 = Message(Univariate, Bernoulli, p=.9)
    msg1 = ruleSPCVIInFactorX(node.id, 1, msg_out, msg_in1, msg_in2, msg_in3)
    msgout = ruleSPCVIOutNFactorNodeX(node.id, nothing, msg_in1, msg_in2, msg_in3)
    @test node.q_memory != node.q
    prod!(msgout.dist, msg_out.dist)
    @test node.q_memory == node.q

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    w = Variable()
    y = Variable()
    node =  Cvi(y, x, z, w, g=f4, opt=[opt,opt,opt], num_samples=1000, num_iterations=[1,1,1], q=[vague(GaussianMeanVariance), vague(Gamma), vague(Bernoulli)], online_inference=[false,true,false])
    msg_out = Message(Univariate, GaussianMeanPrecision, m=1.5, w=0.2)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny)
    msg_in2 = Message(Univariate, Gamma, a=1., b=1.)
    msg_in3 = Message(Univariate, Bernoulli, p=.9)
    msg1 = ruleSPCVIInFactorX(node.id, 1, msg_out, msg_in1, msg_in2, msg_in3)
    msg2 = ruleSPCVIInFactorX(node.id, 2, msg_out, msg_in1, msg_in2, msg_in3)
    msg3 = ruleSPCVIInFactorX(node.id, 3, msg_out, msg_in1, msg_in2, msg_in3)
    msgout = ruleSPCVIOutNFactorNodeX(node.id, nothing, msg_in1, msg_in2, msg_in3)
    @test node.q_memory != node.q
    prod!(msgout.dist, msg_out.dist)
    @test node.q_memory == node.q

    fg = FactorGraph()
    x = Variable()
    z = Variable()
    w = Variable()
    y = Variable()
    node =  Cvi(y, x, z, w, g=f4, opt=[opt,opt,opt], num_samples=1000, num_iterations=[1,1,1], q=[vague(GaussianMeanVariance), vague(Gamma), vague(Bernoulli)], online_inference=[false,true,false])
    msg_out = Message(Univariate, GaussianMeanPrecision, m=1.5, w=0.2)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny)
    msg_in2 = Message(Univariate, Gamma, a=1., b=1.)
    msg_in3 = Message(Univariate, Bernoulli, p=.9)
    msg1 = ruleSPCVIInFactorX(node.id, 1, msg_out, msg_in1, msg_in2, msg_in3)
    msgout = ruleSPCVIOutNFactorNodeX(node.id, nothing, msg_in1, msg_in2, msg_in3)
    @test node.q_memory != node.q
    prod!(msgout.dist, msg_out.dist)
    @test node.q_memory == node.q
    msg2 = ruleSPCVIInFactorX(node.id, 2, msg_out, msg_in1, msg_in2, msg_in3)
    @test node.q_memory == node.q
    msg3 = ruleSPCVIInFactorX(node.id, 3, msg_out, msg_in1, msg_in2, msg_in3)
    @test node.q_memory == node.q

    msg_out = Message(Univariate, Beta, a=1.5, b=0.2)
    @test isUnivariateSampleList(prod!(msgout.dist, msg_out.dist))
end

@testset "MCVIFactorX and JointIndependentProbDist" begin
    fg = FactorGraph()
    x = Variable()
    z = Variable()
    w = Variable()
    y = Variable()
    q1 = ProbabilityDistribution(Univariate,GaussianMeanVariance,m=1.5,v=2.)
    q2 = ProbabilityDistribution(Univariate,Gamma,a=2.5,b=5.)
    q3 = ProbabilityDistribution(Univariate,Bernoulli,p=.45)
    q_list = [q1, q2, q3]
    node =  Cvi(y, x, z, w, g=f4, opt=[opt,opt,opt], num_samples=1000, num_iterations=[1,1,1], q=q_list, online_inference=[false,true,false])
    msg_out = Message(Univariate, GaussianMeanPrecision, m=1.5, w=0.2)
    msg_in1 = Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny)
    msg_in2 = Message(Univariate, Gamma, a=1., b=1.)
    msg_in3 = Message(Univariate, Bernoulli, p=.9)
    joint = ruleMCVIFactorX(node.id, msg_out, msg_in1, msg_in2, msg_in3)
    @test differentialEntropy(joint) == differentialEntropy(q1) + differentialEntropy(q2) + differentialEntropy(q3)
end

end