export
ruleVBWishartOut,
ruleSPWishartOutVPP

ruleSPWishartOutVPP(msg_out::Void, 
                    msg_v::Message{PointMass, MatrixVariate},
                    msg_nu::Message{PointMass, Univariate}) =
    Message(MatrixVariate, Wishart, v=msg_v.dist.params[:m], nu=msg_nu.dist.params[:m])

ruleVBWishartOut(   dist_out::Any,
                    dist_v::ProbabilityDistribution{MatrixVariate},
                    dist_nu::ProbabilityDistribution{Univariate}) =
    Message(MatrixVariate, Wishart, v=unsafeMean(dist_v), nu=unsafeMean(dist_nu))