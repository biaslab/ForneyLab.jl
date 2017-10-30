export  
ruleSPAdditionOutGG,
ruleSPAdditionOutGP,
ruleSPAdditionOutPG,
ruleSPAdditionIn1GG,
ruleSPAdditionIn1PG,
ruleSPAdditionIn2GG

function ruleSPAdditionOutGG(   msg_out::Void,
                                msg_in1::Message{Univariate{Gaussian}},
                                msg_in2::Message{Univariate{Gaussian}})

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v] + msg_in2.dist.params[:v]))
end

function ruleSPAdditionIn2GG(   msg_out::Message{Univariate{Gaussian}},
                                msg_in1::Message{Univariate{Gaussian}},
                                msg_in2::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m], v=msg_out.dist.params[:v] + msg_in1.dist.params[:v]))
end

ruleSPAdditionIn1GG(msg_out::Message{Univariate{Gaussian}}, ::Void, msg_in2::Message{Univariate{Gaussian}}) = ruleSPAdditionIn2GG(msg_out, msg_in2, nothing)

# TODO: add other combinations
function ruleSPAdditionOutGP(   msg_out::Void,
                                msg_in1::Message{Univariate{Gaussian}},
                                msg_in2::Message{Univariate{PointMass}})

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v]))
end
ruleSPAdditionOutPG(::Void, msg_in1::Message{Univariate{PointMass}}, msg_in2::Message{Univariate{Gaussian}}) = ruleSPAdditionOutGP(nothing, msg_in2, msg_in1)

function ruleSPAdditionIn1PG(   msg_out::Message{Univariate{PointMass}},
                                msg_in1::Void,
                                msg_in2::Message{Univariate{Gaussian}})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_in2.dist.params[:v]))
end
