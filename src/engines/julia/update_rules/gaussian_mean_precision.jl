export
ruleVBGaussianMeanPrecisionM, 
ruleVBGaussianMeanPrecisionW, 
ruleVBGaussianMeanPrecisionOut

ruleVBGaussianMeanPrecisionM(   dist_out::ProbabilityDistribution{Univariate},
                                dist_mean::Any,
                                dist_prec::ProbabilityDistribution{Univariate}) =
    Message(Univariate(Gaussian, m=unsafeMean(dist_out), w=unsafeMean(dist_prec)))

ruleVBGaussianMeanPrecisionW(   dist_out::ProbabilityDistribution{Univariate},
                                dist_mean::ProbabilityDistribution{Univariate},
                                dist_prec::Any) =
    Message(Univariate(Gamma, a=1.5, b=0.5*(unsafeVar(dist_mean) + unsafeVar(dist_out) + (unsafeMean(dist_mean) - unsafeMean(dist_out))^2)))

ruleVBGaussianMeanPrecisionOut( dist_out::Any,
                                dist_mean::ProbabilityDistribution{Univariate},
                                dist_prec::ProbabilityDistribution{Univariate}) =
    Message(Univariate(Gaussian, m=unsafeMean(dist_mean), w=unsafeMean(dist_prec)))
