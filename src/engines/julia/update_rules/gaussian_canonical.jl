export
ruleSPGaussianCanonicalOutNPP,
ruleVBGaussianCanonicalOut

ruleSPGaussianCanonicalOutNPP(msg_out::Nothing,
                              msg_weighted_mean::Message{PointMass, V},
                              msg_prec::Message{PointMass}) where V<:VariateType =
    Message(V, Gaussian{Canonical}, xi=deepcopy(msg_weighted_mean.dist.params[:m]), w=deepcopy(msg_prec.dist.params[:m]))

ruleVBGaussianCanonicalOut(dist_out::Any,
                           dist_weighted_mean::Distribution{V},
                           dist_prec::Distribution) where V<:VariateType =
    Message(V, Gaussian{Canonical}, xi=unsafeMean(dist_weighted_mean), w=unsafeMean(dist_prec))