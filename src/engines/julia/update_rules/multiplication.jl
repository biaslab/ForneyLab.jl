export
ruleSPMultiplicationOutVPG,
ruleSPMultiplicationOutVGP,
ruleSPMultiplicationIn1GVP,
ruleSPMultiplicationIn2GPV

function ruleSPMultiplicationOutVPG(msg_out::Void,
                                    msg_in1::Message{PointMass, Univariate},
                                    msg_in2::Message{Gaussian, Univariate})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_in1.dist.params[:m]*msg_in2.dist.params[:m], v=msg_in1.dist.params[:m]^2 * msg_in2.dist.params[:v])
end
ruleSPMultiplicationOutVGP(msg_out::Void, msg_in1::Message{Gaussian, Univariate}, msg_in2::Message{PointMass, Univariate}) = ruleSPMultiplicationOutVPG(nothing, msg_in2, msg_in1)

function ruleSPMultiplicationIn1GVP(msg_out::Message{Gaussian, Univariate},
                                    msg_in1::Void,
                                    msg_in2::Message{PointMass, Univariate})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_out.dist.params[:m]/msg_in2.dist.params[:m], v=msg_out.dist.params[:v]/msg_in2.dist.params[:m]^2)
end
ruleSPMultiplicationIn2GPV(msg_out::Message{Gaussian, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_in2::Void) = ruleSPMultiplicationIn1GVP(msg_out, nothing, msg_in1)