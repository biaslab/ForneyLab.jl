export ruleVBGammaOut

ruleVBGammaOut( dist_out::Any,
                dist_a::ProbabilityDistribution{Univariate},
                dist_b::ProbabilityDistribution{Univariate}) =
    Message(Univariate(Gamma, a=unsafeMean(dist_a), b=unsafeMean(dist_b)))