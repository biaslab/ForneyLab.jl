export
ruleSPDirichletOutVP,
ruleVBDirichletOut

ruleSPDirichletOutVP{V<:VariateType}(msg_out::Void, msg_a::Message{PointMass, V}) = Message(V, Dirichlet, a=deepcopy(msg_a.dist.params[:m]))

ruleVBDirichletOut{V<:VariateType}(marg_out::Any, marg_a::ProbabilityDistribution{V, PointMass}) = Message(V, Dirichlet, a=unsafeMean(marg_a))