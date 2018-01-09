export 
ruleVBLogNormalOut,
ruleSPLogNormalOutVPP

ruleSPLogNormalOutVPP(  msg_out::Void, 
                        msg_m::Message{PointMass, Univariate},
                        msg_s::Message{PointMass, Univariate}) =
    Message(Univariate, LogNormal, m=msg_m.dist.params[:m], s=msg_s.dist.params[:m])

ruleVBLogNormalOut( dist_out::Any,
                    dist_m::ProbabilityDistribution{Univariate},
                    dist_s::ProbabilityDistribution{Univariate}) =
    Message(Univariate, LogNormal, m=unsafeMean(dist_m), s=unsafeMean(dist_s))