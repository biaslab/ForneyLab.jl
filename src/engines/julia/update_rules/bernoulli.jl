export ruleSPBernoulliPV, ruleVBBernoulli2

ruleSPBernoulliPV(msg_in::Message{PointMass}, msg_out::Void) = Message(Bernoulli, p=msg_in.dist.params[:m])

function ruleVBBernoulli2(marg_in::ProbabilityDistribution, marg_out::Any)
    rho_1 = clamp(exp(unsafeLogMean(marg_in)), tiny, huge)
    rho_2 = clamp(exp(unsafeMirroredLogMean(marg_in)), tiny, huge)
    
    Message(Bernoulli, p=rho_1/(rho_1 + rho_2))
end