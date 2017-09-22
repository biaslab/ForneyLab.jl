export ruleVBGamma3

ruleVBGamma3(   dist_a::ProbabilityDistribution,
                dist_b::ProbabilityDistribution,
                dist_out::Any) =
    Message(Gamma, a=unsafeMean(dist_a), b=unsafeMean(dist_b))