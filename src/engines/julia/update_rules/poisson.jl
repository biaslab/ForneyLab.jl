export
ruleSPPoissonOutVP,
ruleSPPoissonLPV,
ruleVBPoissonOut,
ruleVBPoissonL

ruleSPPoissonOutVP(msg_out::Nothing, msg_l::Message{PointMass, Univariate}) =
    Message(Univariate, Poisson, l=msg_l.dist.params[:m])

ruleSPPoissonLPV(msg_out::Message{PointMass, Univariate}, msg_l::Nothing) =
    Message(Univariate, Gamma, a=msg_out.dist.params[:m] + 1.0, b=1.0)

ruleVBPoissonOut(marg_out::Any, marg_l::ProbabilityDistribution{Univariate, Gamma}) =
    Message(Univariate, Poisson, l=(marg_l.params[:a] - 0.5) / marg_l.params[:b])

ruleVBPoissonL(marg_out::ProbabilityDistribution{Univariate}, marg_l::Any) =
    Message(Univariate, Gamma, a=unsafeMean(marg_out) + 1.0, b=1.0)
