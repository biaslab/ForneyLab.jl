export
ruleSPGaussianMeanVarianceOutVPP,
ruleSPGaussianMeanVarianceMPVP,
ruleSPGaussianMeanVarianceOutVGP, 
ruleSPGaussianMeanVarianceMGVP, 
ruleVBGaussianMeanVarianceM,
ruleVBGaussianMeanVarianceOut

ruleSPGaussianMeanVarianceOutVPP{V<:VariateType}(   msg_out::Void,
                                                    msg_mean::Message{PointMass, V},
                                                    msg_var::Message{PointMass}) =
    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), v=deepcopy(msg_var.dist.params[:m]))

ruleSPGaussianMeanVarianceMPVP(msg_out::Message{PointMass}, msg_mean::Void, msg_var::Message{PointMass}) = ruleSPGaussianMeanVarianceOutVPP(msg_mean, msg_out, msg_var)

function ruleSPGaussianMeanVarianceOutVGP{V<:VariateType}(  msg_out::Void,
                                                            msg_mean::Message{Gaussian, V},
                                                            msg_var::Message{PointMass})

    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), v=msg_mean.dist.params[:v] + msg_var.dist.params[:m])
end

ruleSPGaussianMeanVarianceMGVP(msg_out::Message{Gaussian}, msg_mean::Void, msg_var::Message{PointMass}) = ruleSPGaussianMeanVarianceOutVGP(msg_mean, msg_out, msg_var)

ruleVBGaussianMeanVarianceM{V<:VariateType}(dist_out::ProbabilityDistribution{V},
                                            dist_mean::Any,
                                            dist_var::ProbabilityDistribution) =
    Message(V, Gaussian, m=unsafeMean(dist_out), v=unsafeMean(dist_var))

ruleVBGaussianMeanVarianceOut{V<:VariateType}(  dist_out::Any,
                                                dist_mean::ProbabilityDistribution{V},
                                                dist_var::ProbabilityDistribution) =
    Message(V, Gaussian, m=unsafeMean(dist_mean), v=unsafeMean(dist_var))