export
ruleSPGaussianMeanVarianceOutPP,
ruleSPGaussianMeanVarianceMPP,
ruleSPGaussianMeanVarianceOutGP, 
ruleSPGaussianMeanVarianceMGP, 
ruleVBGaussianMeanVarianceM,
ruleVBGaussianMeanVarianceOut

ruleSPGaussianMeanVarianceOutPP(msg_out::Void,
                                msg_mean::Message{PointMass, Univariate},
                                msg_var::Message{PointMass, Univariate}) =
    Message(Univariate(Gaussian, m=msg_mean.dist.params[:m], v=msg_var.dist.params[:m]))

ruleSPGaussianMeanVarianceMPP(msg_out::Message{PointMass, Univariate}, msg_mean::Void, msg_var::Message{PointMass, Univariate}) = ruleSPGaussianMeanVarianceOutPP(msg_mean, msg_out, msg_var)

function ruleSPGaussianMeanVarianceOutGP(   msg_out::Void,
                                            msg_mean::Message{Gaussian, Univariate},
                                            msg_var::Message{PointMass, Univariate})

    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_mean.dist.params[:m], v=msg_mean.dist.params[:v] + msg_var.dist.params[:m]))
end

ruleSPGaussianMeanVarianceMGP(msg_out::Message{Gaussian, Univariate}, msg_mean::Void, msg_var::Message{PointMass, Univariate}) = ruleSPGaussianMeanVarianceOutGP(msg_mean, msg_out, msg_var)

ruleVBGaussianMeanVarianceM(dist_out::ProbabilityDistribution{Univariate},
                            dist_mean::Any,
                            dist_var::ProbabilityDistribution{Univariate}) =
    Message(Univariate(Gaussian, m=unsafeMean(dist_out), v=unsafeMean(dist_var)))

ruleVBGaussianMeanVarianceOut(  dist_out::Any,
                                dist_mean::ProbabilityDistribution{Univariate},
                                dist_var::ProbabilityDistribution{Univariate}) =
    Message(Univariate(Gaussian, m=unsafeMean(dist_mean), v=unsafeMean(dist_var)))
