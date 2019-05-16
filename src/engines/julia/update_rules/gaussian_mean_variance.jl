export
ruleSPGaussianMeanVarianceOutNPP,
ruleSPGaussianMeanVarianceMPNP,
ruleSPGaussianMeanVarianceOutNGP, 
ruleSPGaussianMeanVarianceMGNP, 
ruleVBGaussianMeanVarianceM,
ruleVBGaussianMeanVarianceOut

ruleSPGaussianMeanVarianceOutNPP(   msg_out::Nothing,
                                    msg_mean::Message{PointMass, V},
                                    msg_var::Message{PointMass}) where V<:VariateType =
    Message(V, GaussianMeanVariance, m=deepcopy(msg_mean.dist.params[:m]), v=deepcopy(msg_var.dist.params[:m]))

ruleSPGaussianMeanVarianceMPNP(msg_out::Message{PointMass}, msg_mean::Nothing, msg_var::Message{PointMass}) =
    ruleSPGaussianMeanVarianceOutNPP(msg_mean, msg_out, msg_var)

function ruleSPGaussianMeanVarianceOutNGP(  msg_out::Nothing,
                                            msg_mean::Message{F, V},
                                            msg_var::Message{PointMass}) where {F<:Gaussian, V<:VariateType}

    d_mean = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_mean.dist)

    Message(V, GaussianMeanVariance, m=d_mean.params[:m], v=d_mean.params[:v] + msg_var.dist.params[:m])
end

ruleSPGaussianMeanVarianceMGNP(msg_out::Message{F}, msg_mean::Nothing, msg_var::Message{PointMass}) where F<:Gaussian = 
    ruleSPGaussianMeanVarianceOutNGP(msg_mean, msg_out, msg_var)

ruleVBGaussianMeanVarianceM(dist_out::ProbabilityDistribution{V},
                            dist_mean::Any,
                            dist_var::ProbabilityDistribution) where V<:VariateType =
    Message(V, GaussianMeanVariance, m=unsafeMean(dist_out), v=unsafeMean(dist_var))

ruleVBGaussianMeanVarianceOut(  dist_out::Any,
                                dist_mean::ProbabilityDistribution{V},
                                dist_var::ProbabilityDistribution) where V<:VariateType =
    Message(V, GaussianMeanVariance, m=unsafeMean(dist_mean), v=unsafeMean(dist_var))