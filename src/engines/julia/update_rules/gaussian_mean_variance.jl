export ruleSPGaussianMeanVariancePPV, ruleSPGaussianMeanVarianceVPP, ruleSPGaussianMeanVarianceGPV, ruleSPGaussianMeanVarianceVPG, ruleVBGaussianMeanVariance1, ruleVBGaussianMeanVariance3

# TODO: in-place operation on outbound?
ruleSPGaussianMeanVariancePPV(  msg_mean::Message{PointMass},
                                msg_var::Message{PointMass},
                                msg_out::Void) =
    Message(Gaussian, m=msg_mean.dist.params[:m], v=msg_var.dist.params[:m])

ruleSPGaussianMeanVarianceVPP(msg_mean::Void, msg_var::Message{PointMass}, msg_out::Message{PointMass}) = ruleSPGaussianMeanVariancePPV(msg_out, msg_var, msg_mean)

function ruleSPGaussianMeanVarianceGPV( msg_mean::Message{Gaussian},
                                        msg_var::Message{PointMass},
                                        msg_out::Void)

    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(Gaussian, m=msg_mean.dist.params[:m], v=msg_mean.dist.params[:v] + msg_var.dist.params[:m])
end

ruleSPGaussianMeanVarianceVPG(msg_mean::Void, msg_var::Message{PointMass}, msg_out::Message{Gaussian}) = ruleSPGaussianMeanVarianceGPV(msg_out, msg_var, msg_mean)

ruleVBGaussianMeanVariance1(dist_mean::Any,
                            dist_var::ProbabilityDistribution,
                            dist_out::ProbabilityDistribution) =
    Message(Gaussian, m=unsafeMean(dist_out), v=unsafeMean(dist_var))

ruleVBGaussianMeanVariance3(dist_mean::ProbabilityDistribution,
                            dist_var::ProbabilityDistribution,
                            dist_out::Any) =
    Message(Gaussian, m=unsafeMean(dist_mean), v=unsafeMean(dist_var))
