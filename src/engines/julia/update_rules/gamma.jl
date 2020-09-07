export
ruleVBGammaOut,
ruleSPGammaOutNPP,
ruleVBGammaA,
ruleVBGammaB

ruleSPGammaOutNPP(  msg_out::Nothing,
                    msg_a::Message{PointMass, Univariate},
                    msg_b::Message{PointMass, Univariate}) =
    Message(Univariate, Gamma, a=deepcopy(msg_a.dist.params[:m]), b=deepcopy(msg_b.dist.params[:m]))

ruleVBGammaOut( dist_out::Any,
                dist_a::ProbabilityDistribution{Univariate},
                dist_b::ProbabilityDistribution{Univariate}) =
    Message(Univariate, Gamma, a=unsafeMean(dist_a), b=unsafeMean(dist_b))

ruleVBGammaA( dist_out::ProbabilityDistribution{Univariate},
                dist_a::Any,
                dist_b::ProbabilityDistribution{Univariate}) =
    Message(Univariate, Function, log_pdf = (a)->a*unsafeLogMean(dist_b) + (a-1)*unsafeLogMean(dist_out) - logabsgamma(a)[1])

ruleVBGammaB( dist_out::ProbabilityDistribution{Univariate},
              dist_a::ProbabilityDistribution{Univariate},
              dist_b::Any) =
        Message(Univariate, Gamma, a=unsafeMean(dist_a)+1, b=unsafeMean(dist_out))
