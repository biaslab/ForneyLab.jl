export ruleVBGammaOut

ruleVBGammaOut( dist_out::Any,
                dist_a::Univariate,
                dist_b::Univariate) =
    Message(Univariate(Gamma, a=unsafeMean(dist_a), b=unsafeMean(dist_b)))