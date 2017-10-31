export  
ruleSPAdditionOutGG,
ruleSPAdditionOutGP,
ruleSPAdditionOutPG,
ruleSPAdditionIn1GG,
ruleSPAdditionIn1PG,
ruleSPAdditionIn2GG

function ruleSPAdditionOutGG(   msg_out::Void,
                                msg_in1::Message{Gaussian, Univariate},
                                msg_in2::Message{Gaussian, Univariate})

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v] + msg_in2.dist.params[:v]))
end

function ruleSPAdditionIn2GG(   msg_out::Message{Gaussian, Univariate},
                                msg_in1::Message{Gaussian, Univariate},
                                msg_in2::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m], v=msg_out.dist.params[:v] + msg_in1.dist.params[:v]))
end

ruleSPAdditionIn1GG(msg_out::Message{Gaussian, Univariate}, ::Void, msg_in2::Message{Gaussian, Univariate}) = ruleSPAdditionIn2GG(msg_out, msg_in2, nothing)

# TODO: add other combinations
function ruleSPAdditionOutGP(   msg_out::Void,
                                msg_in1::Message{Gaussian, Univariate},
                                msg_in2::Message{PointMass, Univariate})

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v]))
end
ruleSPAdditionOutPG(::Void, msg_in1::Message{PointMass, Univariate}, msg_in2::Message{Gaussian, Univariate}) = ruleSPAdditionOutGP(nothing, msg_in2, msg_in1)

function ruleSPAdditionIn1PG(   msg_out::Message{PointMass, Univariate},
                                msg_in1::Void,
                                msg_in2::Message{Gaussian, Univariate})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_in2.dist.params[:v]))
end
