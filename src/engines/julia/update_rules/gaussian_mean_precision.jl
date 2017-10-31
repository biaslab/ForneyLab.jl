export
ruleVBGaussianMeanPrecisionM, 
ruleVBGaussianMeanPrecisionW, 
ruleVBGaussianMeanPrecisionOut


#############
# Univariate
#############

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


###############
# Multivariate
###############

ruleVBGaussianMeanPrecisionM(   dist_out::ProbabilityDistribution{Multivariate},
                                dist_mean::Any,
                                dist_prec::ProbabilityDistribution{MatrixVariate}) =
    Message(Multivariate(Gaussian, m=unsafeMean(dist_out), w=unsafeMean(dist_prec)))

ruleVBGaussianMeanPrecisionW(   dist_out::ProbabilityDistribution{Multivariate},
                                dist_mean::ProbabilityDistribution{Multivariate},
                                dist_prec::Any) =
    Message(MatrixVariate(Wishart, v=cholinv( unsafeCov(dist_out) + unsafeCov(dist_mean) + (unsafeMean(dist_out) - unsafeMean(dist_mean))*(unsafeMean(dist_out) - unsafeMean(dist_mean))' ), nu=dims(dist_out) + 2.0) )

ruleVBGaussianMeanPrecisionOut( dist_out::Any,
                                dist_mean::ProbabilityDistribution{Multivariate},
                                dist_prec::ProbabilityDistribution{MatrixVariate}) =
    Message(Multivariate(Gaussian, m=unsafeMean(dist_mean), w=unsafeMean(dist_prec)))
