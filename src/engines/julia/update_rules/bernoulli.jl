export
ruleSPBernoulliOutP,
ruleVBBernoulliOut

ruleSPBernoulliOutP(msg_out::Void, msg_in::Message{Univariate{PointMass}}) = Message(Univariate(Bernoulli, p=msg_in.dist.params[:m]))

function ruleVBBernoulliOut(marg_out::Any, marg_in::Univariate)
    rho_1 = clamp(exp(unsafeLogMean(marg_in)), tiny, huge)
    rho_2 = clamp(exp(unsafeMirroredLogMean(marg_in)), tiny, huge)
    
    Message(Univariate(Bernoulli, p=rho_1/(rho_1 + rho_2)))
end