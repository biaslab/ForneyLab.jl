export ruleSPGaussianMeanVariancePPV, ruleSPGaussianMeanVarianceVPP, ruleSPGaussianMeanVarianceGPV, ruleSPGaussianMeanVarianceVPG

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