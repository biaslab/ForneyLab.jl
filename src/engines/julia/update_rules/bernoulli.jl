export
ruleSPBernoulliOutVP,
ruleVBBernoulliOut,
ruleVBBernoulliIn1

ruleSPBernoulliOutVP(msg_out::Void, msg_p::Message{PointMass, Univariate}) = Message(Univariate, Bernoulli, p=msg_p.dist.params[:m])

function ruleVBBernoulliOut(marg_out::Any, marg_p::ProbabilityDistribution{Univariate})
    rho_1 = clamp(exp(unsafeLogMean(marg_p)), tiny, huge)
    rho_2 = clamp(exp(unsafeMirroredLogMean(marg_p)), tiny, huge)
    
    Message(Univariate, Bernoulli, p=rho_1/(rho_1 + rho_2))
end

ruleVBBernoulliIn1(marg_out::ProbabilityDistribution, marg_p::Any) = Message(Univariate, Beta, a=unsafeMean(marg_out) + 1.0, b= 2.0 - unsafeMean(marg_out))