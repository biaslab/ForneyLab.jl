export
ruleSPDirichletOutNP,
ruleVBDirichletOut,
ruleVBDirichletIn1

ruleSPDirichletOutNP(msg_out::Nothing, msg_a::Message{PointMass, V}) where V<:VariateType = Message(V, Dirichlet, a=deepcopy(msg_a.dist.params[:m]))

ruleVBDirichletOut(marg_out::Any, marg_a::ProbabilityDistribution{V}) where V<:VariateType = Message(V, Dirichlet, a=unsafeMean(marg_a))

function ruleVBDirichletIn1(marg_out::ProbabilityDistribution{Multivariate}, marg_a::Any)
    Message(Multivariate, Function, log_pdf = (a)->sum((a.-1).*unsafeLogMean(marg_out)) - sum(loggamma.(a)) + loggamma(sum(a)))
end
