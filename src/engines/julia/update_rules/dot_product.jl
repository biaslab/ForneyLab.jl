export
ruleSPDotProductOutNPG,
ruleSPDotProductOutNGP,
ruleSPDotProductIn2GPN,
ruleSPDotProductIn1GNP

function ruleSPDotProductOutNPG(msg_out::Nothing,
                                msg_in1::Message{PointMass, Multivariate},
                                msg_in2::Message{F, Multivariate}) where F<:Gaussian

    x = msg_in1.dist.params[:m]
    β = convert(Distribution{Multivariate, Gaussian{Moments}}, msg_in2.dist)
    out_m = x' * β.params[:m]
    out_v = x' * β.params[:v] * x

    return Message(Univariate, Gaussian{Moments}, m=out_m, v=out_v)

end


function ruleSPDotProductOutNGP(msg_out::Nothing,
                                msg_in1::Message{F, Multivariate},
                                msg_in2::Message{PointMass, Multivariate}) where F<:Gaussian

    ruleSPDotProductOutNPG(msg_out, msg_in2, msg_in1)
end


function ruleSPDotProductIn2GPN(msg_out::Message{F, Univariate},
                                msg_in1::Message{PointMass, Multivariate},
                                msg_in2::Nothing) where F<:Gaussian

    x = msg_in1.dist.params[:m] # We'll call in1 x
    d = length(x)

    y = convert(Distribution{Univariate, Gaussian{Canonical}}, msg_out.dist)
    xi = x * y.params[:xi]
    w = x * y.params[:w] * x'
    w += tiny*diageye(size(w)[1]) # Ensure w is invertible

    return Message(Multivariate, Gaussian{Canonical}, xi=xi, w=w)
end


function ruleSPDotProductIn1GNP(msg_out::Message{F, Univariate},
                                msg_in1::Nothing,
                                msg_in2::Message{PointMass, Multivariate}) where F<:Gaussian

    ruleSPDotProductIn2GPN(msg_out, msg_in2, msg_in1)
end
