export
ruleSPDirichletOutVP,
ruleVBDirichletOut

ruleSPDirichletOutVP(msg_out::Nothing, msg_a::Message{PointMass, V}) where V<:VariateType = Message(V, Dirichlet, a=deepcopy(msg_a.dist.params[:m]))

ruleVBDirichletOut(marg_out::Any, marg_a::ProbabilityDistribution{V, PointMass}) where V<:VariateType = Message(V, Dirichlet, a=unsafeMean(marg_a))