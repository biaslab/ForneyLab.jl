export ruleVBGaussianMeanPrecision1, ruleVBGaussianMeanPrecision2, ruleVBGaussianMeanPrecision3

ruleVBGaussianMeanPrecision1(   dist_mean::Any,
                                dist_prec::ProbabilityDistribution,
                                dist_out::ProbabilityDistribution) =
    Message(Gaussian, m=unsafeMean(dist_out), w=unsafeMean(dist_prec))

ruleVBGaussianMeanPrecision2(   dist_mean::ProbabilityDistribution,
                                dist_prec::Any,
                                dist_out::ProbabilityDistribution) =
    Message(Gamma, a=1.5, b=0.5*(unsafeVar(dist_mean) + unsafeVar(dist_out) + (unsafeMean(dist_mean) - unsafeMean(dist_out))^2))

ruleVBGaussianMeanPrecision3(   dist_mean::ProbabilityDistribution,
                                dist_prec::ProbabilityDistribution,
                                dist_out::Any) =
    Message(Gaussian, m=unsafeMean(dist_mean), w=unsafeMean(dist_prec))
