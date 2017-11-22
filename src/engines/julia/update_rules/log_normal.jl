export ruleVBLogNormalOut

ruleVBLogNormalOut( dist_out::Any,
                	dist_m::ProbabilityDistribution{Univariate},
                	dist_s::ProbabilityDistribution{Univariate}) =
    Message(Univariate, LogNormal, m=unsafeMean(dist_m), s=unsafeMean(dist_s))