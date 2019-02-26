export
ruleSPBernoulliOutVP,
ruleSPBernoulliIn1PV,
ruleSPBernoulliOutVB,
ruleVBBernoulliOut,
ruleVBBernoulliIn1

ruleSPBernoulliOutVP(msg_out::Nothing, msg_p::Message{PointMass, Univariate}) = Message(Univariate, Bernoulli, p=deepcopy(msg_p.dist.params[:m]))

ruleSPBernoulliIn1PV(msg_out::Message{PointMass, Univariate}, msg_p::Nothing) = Message(Univariate, Beta, a=msg_out.dist.params[:m]+1, b=2-msg_out.dist.params[:m])

ruleSPBernoulliOutVB(msg_out::Nothing, msg_p::Message{Beta, Univariate}) = Message(Univariate, Bernoulli, p=msg_p.dist.params[:a]/(msg_p.dist.params[:a] + msg_p.dist.params[:b]))

function ruleVBBernoulliOut(marg_out::Any, marg_p::ProbabilityDistribution{Univariate})
    rho_1 = clamp(exp(unsafeLogMean(marg_p)), tiny, huge)
    rho_2 = clamp(exp(unsafeMirroredLogMean(marg_p)), tiny, huge)

    Message(Univariate, Bernoulli, p=rho_1/(rho_1 + rho_2))
end

ruleVBBernoulliIn1(marg_out::ProbabilityDistribution, marg_p::Any) = Message(Univariate, Beta, a=unsafeMean(marg_out) + 1.0, b= 2.0 - unsafeMean(marg_out))
