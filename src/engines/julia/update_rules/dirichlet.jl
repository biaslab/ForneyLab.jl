export
ruleSPDirichletOutVP,
ruleVBDirichletOut

ruleSPDirichletOutVP(msg_out::Void, msg_a::Message{PointMass, Multivariate}) = Message(Multivariate, Dirichlet, a=deepcopy(msg_a.dist.params[:m]))

ruleVBDirichletOut(marg_out::Any, marg_a::ProbabilityDistribution{Multivariate, PointMass}) = Message(Multivariate, Dirichlet, a=unsafeMean(marg_a))