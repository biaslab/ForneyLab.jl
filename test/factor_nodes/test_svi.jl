module SviTest
    
using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeVar
using ForneyLab: SPSVIOutNM, SPSVIIn1MN
using Flux.Optimise

@testset "Svi Node Construction" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = ForgetDelayDescent(1.,0.7)
    @test Svi(y, x, opt=opt, q=vague(GaussianMeanVariance), batch_size=1, dataset_size=10).q_memory == vague(GaussianMeanVariance)
end

#-------------
# Update rules
#-------------

# SVI has to be used within VMP so the naming with SP just refers to the fact that we designed an SP like updates around this node
@testset "SPSVIOutNM" begin
    @test SPSVIOutNM <: SumProductRule{SVI}
    @test outboundType(SPSVIOutNM) == Message{SetProbDist}
    @test isApplicable(SPSVIOutNM, [Nothing, Message{Gaussian}])
    @test isApplicable(SPSVIOutNM, [Nothing, Message{Bernoulli}])

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = ForgetDelayDescent(1.,0.7)
    node = Svi(y, x, opt=opt, q=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2., v=4.), batch_size=1, dataset_size=10)
    msg = ruleSPSVIOutNM(node.id,nothing, Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny))

    @test msg.dist.params[:message] == Message(Univariate, GaussianMeanVariance, m=2.0, v=tiny)
    @test msg.dist.params[:q] == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2., v=4.)
    @test msg.dist.params[:node_id] == node.id
end

@testset "SPSVIIn1MN" begin
    @test SPSVIIn1MN <: SumProductRule{SVI}
    @test outboundType(SPSVIIn1MN) == Message{FactorNode}
    @test isApplicable(SPSVIIn1MN, [Message{Gaussian}, Nothing])
    @test isApplicable(SPSVIIn1MN, [Message{Gaussian}, Nothing])

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2., v=4.), batch_size=1, dataset_size=1)
    msg_in = Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    msg_out = Message(Univariate, GaussianMeanPrecision, m=2.0, w=2.0)
    msg = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    @test isapprox(unsafeMean(msg.dist), unsafeMean(msg_out.dist))
    @test isapprox(unsafeVar(msg.dist), unsafeVar(msg_out.dist))

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2., v=4.), batch_size=1, dataset_size=10)
    msg_in = Message(Univariate, GaussianMeanVariance, m=1.0, v=3.0)
    msg_out = Message(Univariate, GaussianMeanPrecision, m=2.0, w=2.0)
    msg = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    prod_res = msg_out.dist
    for _=1:9 prod_res = prod!(prod_res,msg_out.dist) end
    @test prod_res == msg.dist

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=vague(Gamma), batch_size=1, dataset_size=1)
    msg_in = Message(Univariate, Gamma, a=1.0, b=3.0)
    msg_out = Message(Univariate, Gamma, a=2.0, b=2.0)
    msg = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    @test msg == msg_out

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=vague(Wishart,2), batch_size=1, dataset_size=1)
    msg_in = Message(MatrixVariate, Wishart, v=[1.0 0.;0. 1.], nu=2.0)
    msg_out = Message(MatrixVariate, Wishart, v=[1.5 0.2;0.2 1.], nu=4.0)
    msg = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    @test msg == msg_out
end

#-------------
# SetProbDist
#-------------

@testset "SetProbDist" begin
    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=2., w=4.), batch_size=1, dataset_size=1)
    @test node.q == node.q_memory
    msg_in = Message(Univariate, GaussianMeanPrecision, m=1.0, w=3.0)
    msg_out = Message(Univariate, GaussianWeightedMeanPrecision, xi=2.0, w=2.0)
    # msg_out = Message(Univariate, GaussianMeanPrecision, m=2.0, w=2.0) # returns ambiguity error ?
    msg_backward = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    @test node.q != node.q_memory
    msg_forward = ruleSPSVIOutNM(node.id, nothing, msg_in)
    prod!(msg_forward.dist, msg_out.dist)
    @test node.q == node.q_memory

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=vague(GaussianMeanVariance,2), batch_size=1, dataset_size=1)
    @test node.q == node.q_memory
    msg_in = Message(Multivariate, GaussianMeanVariance, m=[1.0,2.5], v=[1. 0.;0. 1.])
    msg_out = Message(Multivariate, GaussianMeanVariance, m=[1.0,2.5], v=[1. 0.;0. 1.])
    msg_backward = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    @test node.q != node.q_memory
    msg_forward = ruleSPSVIOutNM(node.id, nothing, msg_in)
    prod!(msg_forward.dist, msg_out.dist)
    @test node.q == node.q_memory

    fg = FactorGraph()
    x = Variable()
    y = Variable()
    opt = Descent(1.0)
    node = Svi(y, x, opt=opt, q=vague(Gamma), batch_size=1, dataset_size=1)
    @test node.q == node.q_memory
    msg_in = Message(Univariate, Gamma, a=1.0, b=3.0)
    msg_out = Message(Univariate, Gamma, a=2.0, b=2.0)
    msg_backward = ruleSPSVIIn1MN(node.id, msg_out, msg_in)
    @test node.q != node.q_memory
    msg_forward = ruleSPSVIOutNM(node.id, nothing, msg_in)
    prod!(msg_forward.dist, msg_out.dist)
    @test node.q == node.q_memory
end

end