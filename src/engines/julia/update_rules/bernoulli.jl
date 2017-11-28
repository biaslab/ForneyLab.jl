export
ruleSPBernoulliOutVP,
ruleVBBernoulliOut

ruleSPBernoulliOutVP(msg_out::Void, msg_p::Message{PointMass, Univariate}) = Message(Univariate, Bernoulli, p=msg_p.dist.params[:m])

function ruleVBBernoulliOut(marg_out::Any, marg_p::ProbabilityDistribution{Univariate})
    rho_1 = clamp(exp(unsafeLogMean(marg_p)), tiny, huge)
    rho_2 = clamp(exp(unsafeMirroredLogMean(marg_p)), tiny, huge)
    
    Message(Univariate, Bernoulli, p=rho_1/(rho_1 + rho_2))
end