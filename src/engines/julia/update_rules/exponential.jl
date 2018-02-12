export
ruleSPExponentialOutVG,
ruleSPExponentialOutVP,
ruleSPExponentialIn1LV,
ruleSPExponentialIn1PV

function ruleSPExponentialOutVG(msg_out::Void, 
                                msg_in1::Message{Gaussian, Univariate})

    ensureParameters!(msg_in1.dist, (:m, :v))

    return Message(Univariate, LogNormal, m=deepcopy(msg_in1.dist.params[:m]), s=deepcopy(msg_in1.dist.params[:v]))
end

ruleSPExponentialIn1LV(msg_out::Message{LogNormal, Univariate}, msg_in1::Void) = Message(Univariate, Gaussian, m=msg_out.dist.params[:m], v=msg_out.dist.params[:s])

ruleSPExponentialOutVP(msg_out::Void, msg_in1::Message{PointMass, Univariate}) = Message(Univariate, PointMass, m=exp(msg_in1.dist.params[:m]))

ruleSPExponentialIn1PV(msg_out::Message{PointMass, Univariate}, msg_in1::Void) = Message(Univariate, PointMass, m=log(msg_out.dist.params[:m]))