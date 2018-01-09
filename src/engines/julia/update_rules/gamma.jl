export 
ruleVBGammaOut,
ruleSPGammaOutVPP

ruleSPGammaOutVPP(  msg_out::Void, 
                    msg_a::Message{PointMass, Univariate},
                    msg_b::Message{PointMass, Univariate}) =
    Message(Univariate, Gamma, a=msg_a.dist.params[:m], b=msg_b.dist.params[:m])

ruleVBGammaOut( dist_out::Any,
                dist_a::ProbabilityDistribution{Univariate},
                dist_b::ProbabilityDistribution{Univariate}) =
    Message(Univariate, Gamma, a=unsafeMean(dist_a), b=unsafeMean(dist_b))