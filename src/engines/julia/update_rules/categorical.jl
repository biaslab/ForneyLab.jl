export
ruleSPCategoricalOutNP,
ruleSPCategoricalIn1PN,
ruleVBCategoricalOut,
ruleVBCategoricalIn1

ruleSPCategoricalOutNP(msg_out::Nothing, msg_p::Message{PointMass, Multivariate}) = Message(Univariate, Categorical, p=deepcopy(msg_p.dist.params[:m]))

ruleSPCategoricalIn1PN(msg_out::Message{PointMass, Multivariate}, msg_p::Nothing) = Message(Multivariate, Dirichlet, a=deepcopy(msg_out.dist.params[:m]).+ 1.0)

function ruleVBCategoricalOut(marg_out::Any, marg_p::ProbabilityDistribution{Multivariate})
    rho = clamp.(exp.(unsafeLogMean(marg_p)), tiny, Inf) # Softens the parameter
    
    Message(Univariate, Categorical, p=rho./sum(rho))
end

ruleVBCategoricalIn1(marg_out::ProbabilityDistribution, marg_p::Any) = Message(Multivariate, Dirichlet, a=unsafeMeanVector(marg_out) .+ 1.0)
