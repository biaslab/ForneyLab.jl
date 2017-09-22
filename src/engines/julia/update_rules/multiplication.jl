export ruleSPMultiplicationGPV, ruleSPMultiplicationVPG

function ruleSPMultiplicationGPV(   msg_in::Message{Gaussian},
                                    msg_gain::Message{PointMass},
                                    msg_out::Void)

    ensureParameters!(msg_in.dist, (:m, :v))

    Message(Gaussian, m=msg_gain.dist.params[:m]*msg_in.dist.params[:m], v=msg_gain.dist.params[:m]^2 * msg_in.dist.params[:v])
end

function ruleSPMultiplicationVPG(   msg_in::Void,
                                    msg_gain::Message{PointMass},
                                    msg_out::Message{Gaussian})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Gaussian, m=msg_out.dist.params[:m]/msg_gain.dist.params[:m], v=msg_out.dist.params[:v]/msg_gain.dist.params[:m]^2)
end