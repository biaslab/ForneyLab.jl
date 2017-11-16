export  
ruleSPAdditionOutVGG,
ruleSPAdditionOutVGP,
ruleSPAdditionOutVPG,
ruleSPAdditionIn1GVG,
ruleSPAdditionIn1PVG,
ruleSPAdditionIn1GVP,
ruleSPAdditionIn2GGV,
ruleSPAdditionIn2PGV,
ruleSPAdditionIn2GPV

function ruleSPAdditionOutVGG(   msg_out::Void,
                                msg_in1::Message{Gaussian, Univariate},
                                msg_in2::Message{Gaussian, Univariate})

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v] + msg_in2.dist.params[:v])
end

function ruleSPAdditionIn2GGV(   msg_out::Message{Gaussian, Univariate},
                                msg_in1::Message{Gaussian, Univariate},
                                msg_in2::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m], v=msg_out.dist.params[:v] + msg_in1.dist.params[:v])
end

ruleSPAdditionIn1GVG(msg_out::Message{Gaussian, Univariate}, ::Void, msg_in2::Message{Gaussian, Univariate}) = ruleSPAdditionIn2GGV(msg_out, msg_in2, nothing)

# TODO: add other combinations
function ruleSPAdditionOutVGP(   msg_out::Void,
                                msg_in1::Message{Gaussian, Univariate},
                                msg_in2::Message{PointMass, Univariate})

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v])
end
ruleSPAdditionOutVPG(::Void, msg_in1::Message{PointMass, Univariate}, msg_in2::Message{Gaussian, Univariate}) = ruleSPAdditionOutVGP(nothing, msg_in2, msg_in1)

function ruleSPAdditionIn1PVG(  msg_out::Message{PointMass, Univariate},
                                msg_in1::Void,
                                msg_in2::Message{Gaussian, Univariate})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_in2.dist.params[:v])
end
ruleSPAdditionIn2PGV(msg_out::Message{PointMass, Univariate}, msg_in1::Message{Gaussian, Univariate}, msg_in2::Void) = ruleSPAdditionIn1PVG(msg_out, nothing, msg_in1)

function ruleSPAdditionIn1GVP(  msg_out::Message{Gaussian, Univariate},
                                msg_in1::Void,
                                msg_in2::Message{PointMass, Univariate})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_out.dist.params[:v])
end
ruleSPAdditionIn2GPV(msg_out::Message{Gaussian, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_in2::Void) = ruleSPAdditionIn1GVP(msg_out, nothing, msg_in1)