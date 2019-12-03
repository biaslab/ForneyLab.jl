export
ruleSPGaussianWeightedMeanPrecisionOutNPP,
ruleVBGaussianWeightedMeanPrecisionOut

ruleSPGaussianWeightedMeanPrecisionOutNPP(  msg_out::Nothing,
                                            msg_weighted_mean::Message{PointMass, V},
                                            msg_prec::Message{PointMass}) where V<:VariateType =
    Message(V, GaussianWeightedMeanPrecision, xi=deepcopy(msg_weighted_mean.dist.params[:m]), w=deepcopy(msg_prec.dist.params[:m]))

ruleVBGaussianWeightedMeanPrecisionOut( dist_out::Any,
                                        dist_weighted_mean::ProbabilityDistribution{V},
                                        dist_prec::ProbabilityDistribution) where V<:VariateType =
    Message(V, GaussianWeightedMeanPrecision, xi=unsafeMean(dist_weighted_mean), w=unsafeMean(dist_prec))