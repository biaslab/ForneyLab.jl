export
ruleSPDotProductOutVPG,
ruleSPDotProductOutVGP,
ruleSPDotProductIn2GPV,
ruleSPDotProductIn1GVP

function ruleSPDotProductOutVPG{F<:Gaussian}(   msg_out::Nothing,
                                                msg_in1::Message{PointMass, Multivariate},
                                                msg_in2::Message{F, Multivariate})

    x = msg_in1.dist.params[:m]
    β = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in2.dist)
    out_m = x' * β.params[:m]
    out_v = x' * β.params[:v] * x

    return Message(Univariate, GaussianMeanVariance, m=out_m, v=out_v)

end


function ruleSPDotProductOutVGP{F<:Gaussian}(   msg_out::Nothing,
                                                msg_in1::Message{F, Multivariate},
                                                msg_in2::Message{PointMass, Multivariate})

    ruleSPDotProductOutVPG(msg_out, msg_in2, msg_in1)
end


function ruleSPDotProductIn2GPV{F<:Gaussian}(   msg_out::Message{F, Univariate},
                                                msg_in1::Message{PointMass, Multivariate},
                                                msg_in2::Nothing)

    x = msg_in1.dist.params[:m] # We'll call in1 x
    d = length(x)

    y = convert(ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, msg_out.dist)
    xi = x * y.params[:xi]
    w = x * y.params[:w] * x'
    w += tiny*diageye(size(w)[1]) # Ensure w is invertible

    return Message(Multivariate, GaussianWeightedMeanPrecision, xi=xi, w=w)
end


function ruleSPDotProductIn1GVP{F<:Gaussian}(   msg_out::Message{F, Univariate},
                                                msg_in1::Nothing,
                                                msg_in2::Message{PointMass, Multivariate})

    ruleSPDotProductIn2GPV(msg_out, msg_in2, msg_in1)
end
