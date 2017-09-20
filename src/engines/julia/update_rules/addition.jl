export ruleSPAdditionGGV

function ruleSPAdditionGGV( msg_in1::Message{Gaussian},
                            msg_in2::Message{Gaussian},
                            msg_out::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v] + msg_in2.dist.params[:v])
end