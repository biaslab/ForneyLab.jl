export
ruleSPBetaOutVPP,
ruleVBBetaOut

ruleSPBetaOutVPP(msg_out::Nothing, msg_a::Message{PointMass, Univariate}, msg_b::Message{PointMass, Univariate}) = Message(Univariate, Beta, a=deepcopy(msg_a.dist.params[:m]), b=deepcopy(msg_b.dist.params[:m]))

ruleVBBetaOut(marg_out::Any, marg_a::ProbabilityDistribution{Univariate, PointMass}, marg_b::ProbabilityDistribution{Univariate, PointMass}) = Message(Univariate, Beta, a=marg_a.params[:m], b=marg_b.params[:m])