export
ruleSPDotProductOutVPG,
ruleSPDotProductOutVGP,
ruleSPDotProductIn2GPV,
ruleSPDotProductIn1GVP

function ruleSPDotProductOutVPG(msg_out::Void,
                                msg_in1::Message{PointMass, Multivariate},
                                msg_in2::Message{Gaussian, Multivariate})

    x = msg_in1.dist.params[:m]
    β = ensureParameters!(msg_in2.dist, (:m, :v))
    out_m = x' * β.params[:m]
    out_v = x' * β.params[:v] * x

    return Message(Univariate, Gaussian, m=out_m, v=out_v)

end


function ruleSPDotProductOutVGP(msg_out::Void,
                                msg_in1::Message{Gaussian, Multivariate},
                                msg_in2::Message{PointMass, Multivariate})

    ruleSPDotProductOutVPG(msg_out, msg_in2, msg_in1)
end


function ruleSPDotProductIn2GPV(msg_out::Message{Gaussian, Univariate},
                                msg_in1::Message{PointMass, Multivariate},
                                msg_in2::Void)

    x = msg_in1.dist.params[:m] # We'll call in1 x
    d = length(x)

    y = ensureParameters!(msg_out.dist, (:xi, :w))
    xi = x * y.params[:xi]
    w = x * y.params[:w] * x'

    return Message(Multivariate, Gaussian, xi=xi, w=w)
end


function ruleSPDotProductIn1GVP(msg_out::Message{Gaussian, Univariate},
                                msg_in1::Void,
                                msg_in2::Message{PointMass, Multivariate})

    ruleSPDotProductIn2GPV(msg_out, msg_in2, msg_in1)
end
