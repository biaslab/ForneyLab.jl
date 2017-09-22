export ruleSPAdditionGGV, ruleSPAdditionGVG, ruleSPAdditionVGG, ruleSPAdditionGPV, ruleSPAdditionVGP, ruleSPAdditionPGV

function ruleSPAdditionGGV( msg_in1::Message{Gaussian},
                            msg_in2::Message{Gaussian},
                            msg_out::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v] + msg_in2.dist.params[:v])
end

function ruleSPAdditionGVG( msg_in1::Message{Gaussian},
                            msg_in2::Void,
                            msg_out::Message{Gaussian})

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Gaussian, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m], v=msg_out.dist.params[:v] + msg_in1.dist.params[:v])
end

ruleSPAdditionVGG(msg_in1::Void, msg_in2::Message{Gaussian}, msg_out::Message{Gaussian}) = ruleSPAdditionGVG(msg_in2, msg_in1, msg_out)

# TODO: add other combinations
function ruleSPAdditionGPV( msg_in1::Message{Gaussian},
                            msg_in2::Message{PointMass},
                            msg_out::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v])
end
ruleSPAdditionPGV(msg_in1::Message{PointMass}, msg_in2::Message{Gaussian}, msg_out::Void) = ruleSPAdditionGPV(msg_in2, msg_in1, msg_out)

function ruleSPAdditionVGP( msg_in1::Void,
                            msg_in2::Message{Gaussian},
                            msg_out::Message{PointMass})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_in2.dist.params[:v])
end
