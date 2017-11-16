export
ruleSPCategoricalOutVP,
ruleVBCategoricalOut

ruleSPCategoricalOutVP(msg_out::Void, msg_p::Message{PointMass, Multivariate}) = Message(Univariate, Categorical, p=msg_p.dist.params[:m])

function ruleVBCategoricalOut(marg_out::Any, marg_p::ProbabilityDistribution{Multivariate})
    rho = clamp(exp(unsafeLogMean(marg_p)), tiny, huge)
    
    Message(Univariate, Categorical, p=rho./sum(rho))
end