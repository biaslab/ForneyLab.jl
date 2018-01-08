export
ruleSPGaussianMeanPrecisionOutVPP,
ruleSPGaussianMeanPrecisionMPVP,
ruleSPGaussianMeanPrecisionOutVGP, 
ruleSPGaussianMeanPrecisionMGVP, 
ruleVBGaussianMeanPrecisionM, 
ruleVBGaussianMeanPrecisionW, 
ruleVBGaussianMeanPrecisionOut

ruleSPGaussianMeanPrecisionOutVPP{V<:VariateType}(  msg_out::Void,
                                                    msg_mean::Message{PointMass, V},
                                                    msg_prec::Message{PointMass}) =
    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), w=deepcopy(msg_prec.dist.params[:m]))

ruleSPGaussianMeanPrecisionMPVP(msg_out::Message{PointMass}, msg_mean::Void, msg_prec::Message{PointMass}) = ruleSPGaussianMeanPrecisionOutVPP(msg_mean, msg_out, msg_prec)

function ruleSPGaussianMeanPrecisionOutVGP{V<:VariateType}( msg_out::Void,
                                                            msg_mean::Message{Gaussian, V},
                                                            msg_prec::Message{PointMass})

    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), v=msg_mean.dist.params[:v] + cholinv(msg_prec.dist.params[:m]))
end

ruleSPGaussianMeanPrecisionMGVP(msg_out::Message{Gaussian}, msg_mean::Void, msg_prec::Message{PointMass}) = ruleSPGaussianMeanPrecisionOutVGP(msg_mean, msg_out, msg_prec)

ruleVBGaussianMeanPrecisionM{V<:VariateType}(   dist_out::ProbabilityDistribution{V},
                                                dist_mean::Any,
                                                dist_prec::ProbabilityDistribution) =
    Message(V, Gaussian, m=unsafeMean(dist_out), w=unsafeMean(dist_prec))

ruleVBGaussianMeanPrecisionW(   dist_out::ProbabilityDistribution{Univariate},
                                dist_mean::ProbabilityDistribution{Univariate},
                                dist_prec::Any) =
    Message(Univariate, Gamma, a=1.5, b=0.5*(unsafeVar(dist_mean) + unsafeVar(dist_out) + (unsafeMean(dist_mean) - unsafeMean(dist_out))^2))

ruleVBGaussianMeanPrecisionW(   dist_out::ProbabilityDistribution{Multivariate},
                                dist_mean::ProbabilityDistribution{Multivariate},
                                dist_prec::Any) =
    Message(MatrixVariate, Wishart, v=cholinv( unsafeCov(dist_out) + unsafeCov(dist_mean) + (unsafeMean(dist_out) - unsafeMean(dist_mean))*(unsafeMean(dist_out) - unsafeMean(dist_mean))' ), nu=dims(dist_out) + 2.0) 

ruleVBGaussianMeanPrecisionOut{V<:VariateType}( dist_out::Any,
                                                dist_mean::ProbabilityDistribution{V},
                                                dist_prec::ProbabilityDistribution) =
    Message(V, Gaussian, m=unsafeMean(dist_mean), w=unsafeMean(dist_prec))