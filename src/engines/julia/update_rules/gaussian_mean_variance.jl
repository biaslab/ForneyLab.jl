export ruleSPGaussianMeanVariancePPV, ruleSPGaussianMeanVarianceVPP

# TODO: in-place operation on outbound?
ruleSPGaussianMeanVariancePPV(  msg_mean::Message{PointMass},
                                msg_var::Message{PointMass},
                                msg_out::Void) =
    Message(Gaussian, m=msg_mean.dist.params[:m], v=msg_var.dist.params[:m])

ruleSPGaussianMeanVarianceVPP(  msg_mean::Void,
                                msg_var::Message{PointMass},
                                msg_out::Message{PointMass}) =
    Message(Gaussian, m=msg_out.dist.params[:m], v=msg_var.dist.params[:m])