export  
ruleSPAdditionOutVGG,
ruleSPAdditionOutVGP,
ruleSPAdditionOutVPG,
ruleSPAdditionOutVPP,
ruleSPAdditionIn1GVG,
ruleSPAdditionIn1PVG,
ruleSPAdditionIn1GVP,
ruleSPAdditionIn1PVP,
ruleSPAdditionIn2GGV,
ruleSPAdditionIn2PGV,
ruleSPAdditionIn2GPV,
ruleSPAdditionIn2PPV

function ruleSPAdditionOutVGG{V<:Union{Univariate, Multivariate}}(  msg_out::Void,
                                                                    msg_in1::Message{Gaussian, V},
                                                                    msg_in2::Message{Gaussian, V})

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(V, Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v] + msg_in2.dist.params[:v])
end

function ruleSPAdditionIn2GGV{V<:Union{Univariate, Multivariate}}(  msg_out::Message{Gaussian, V},
                                                                    msg_in1::Message{Gaussian, V},
                                                                    msg_in2::Void)

    ensureParameters!(msg_in1.dist, (:m, :v))
    ensureParameters!(msg_out.dist, (:m, :v))

    Message(V, Gaussian, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m], v=msg_out.dist.params[:v] + msg_in1.dist.params[:v])
end

ruleSPAdditionIn1GVG{V<:Union{Univariate, Multivariate}}(msg_out::Message{Gaussian, V}, ::Void, msg_in2::Message{Gaussian, V}) = ruleSPAdditionIn2GGV(msg_out, msg_in2, nothing)

function ruleSPAdditionOutVGP{V<:Union{Univariate, Multivariate}}(  msg_out::Void,
                                                                    msg_in1::Message{Gaussian, V},
                                                                    msg_in2::Message{PointMass, V})

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(V, Gaussian, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m], v=msg_in1.dist.params[:v])
end
ruleSPAdditionOutVPG{V<:Union{Univariate, Multivariate}}(::Void, msg_in1::Message{PointMass, V}, msg_in2::Message{Gaussian, V}) = ruleSPAdditionOutVGP(nothing, msg_in2, msg_in1)

function ruleSPAdditionIn1PVG{V<:Union{Univariate, Multivariate}}(  msg_out::Message{PointMass, V},
                                                                    msg_in1::Void,
                                                                    msg_in2::Message{Gaussian, V})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(V, Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_in2.dist.params[:v])
end
ruleSPAdditionIn2PGV{V<:Union{Univariate, Multivariate}}(msg_out::Message{PointMass, V}, msg_in1::Message{Gaussian, V}, msg_in2::Void) = ruleSPAdditionIn1PVG(msg_out, nothing, msg_in1)

function ruleSPAdditionIn1GVP{V<:Union{Univariate, Multivariate}}(  msg_out::Message{Gaussian, V},
                                                                    msg_in1::Void,
                                                                    msg_in2::Message{PointMass, V})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(V, Gaussian, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m], v=msg_out.dist.params[:v])
end
ruleSPAdditionIn2GPV{V<:Union{Univariate, Multivariate}}(msg_out::Message{Gaussian, V}, msg_in1::Message{PointMass, V}, msg_in2::Void) = ruleSPAdditionIn1GVP(msg_out, nothing, msg_in1)

ruleSPAdditionOutVPP{V<:Union{Univariate, Multivariate}}(msg_out::Void, msg_in1::Message{PointMass, V}, msg_in2::Message{PointMass, V}) = Message(V, PointMass, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m])
ruleSPAdditionIn2PPV{V<:Union{Univariate, Multivariate}}(msg_out::Message{PointMass, V}, msg_in1::Message{PointMass, V}, msg_in2::Void) = Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m])
ruleSPAdditionIn1PVP{V<:Union{Univariate, Multivariate}}(msg_out::Message{PointMass, V}, msg_in1::Void, msg_in2::Message{PointMass, V}) = Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m])
