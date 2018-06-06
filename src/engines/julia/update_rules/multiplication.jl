export
ruleSPMultiplicationOutVPG,
ruleSPMultiplicationOutVGP,
ruleSPMultiplicationOutVPP,
ruleSPMultiplicationIn1GVP,
ruleSPMultiplicationIn1PVP,
ruleSPMultiplicationAGPV,
ruleSPMultiplicationAPPV

# Univariate (in1 and a commute)
function ruleSPMultiplicationOutVPG(msg_out::Void,
                                    msg_in1::Message{PointMass, Univariate},
                                    msg_a::Message{Gaussian, Univariate})

    ensureParameters!(msg_a.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_in1.dist.params[:m]*msg_a.dist.params[:m], v=msg_a.dist.params[:v]*msg_in1.dist.params[:m]^2)
end
ruleSPMultiplicationOutVGP(msg_out::Void, msg_in1::Message{Gaussian, Univariate}, msg_a::Message{PointMass, Univariate}) = ruleSPMultiplicationOutVPG(nothing, msg_a, msg_in1)

function ruleSPMultiplicationIn1GVP(msg_out::Message{Gaussian, Univariate},
                                    msg_in1::Void,
                                    msg_a::Message{PointMass, Univariate})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_out.dist.params[:m]/msg_a.dist.params[:m], v=msg_out.dist.params[:v]/msg_a.dist.params[:m]^2)
end
ruleSPMultiplicationAGPV(msg_out::Message{Gaussian, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_a::Void) = ruleSPMultiplicationIn1GVP(msg_out, nothing, msg_in1)

ruleSPMultiplicationOutVPP(msg_out::Void, msg_in1::Message{PointMass, Univariate}, msg_a::Message{PointMass, Univariate}) = Message(Univariate, PointMass, m=msg_in1.dist.params[:m]*msg_a.dist.params[:m])
ruleSPMultiplicationIn1PVP(msg_out::Message{PointMass, Univariate}, msg_in1::Void, msg_a::Message{PointMass, Univariate}) = Message(Univariate, PointMass, m=msg_out.dist.params[:m]/msg_a.dist.params[:m])
ruleSPMultiplicationAPPV(msg_out::Message{PointMass, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_a::Void) = Message(Univariate, PointMass, m=msg_out.dist.params[:m]/msg_in1.dist.params[:m])

# Univariate*multivariate (in1 and a commute)
function ruleSPMultiplicationOutVPG(msg_out::Void,
                                    msg_in1::Message{PointMass, Multivariate},
                                    msg_a::Message{Gaussian, Univariate})

    ensureParameters!(msg_a.dist, (:m, :v))

    Message(Multivariate, Gaussian, m=msg_in1.dist.params[:m]*msg_a.dist.params[:m], v=msg_a.dist.params[:v]*msg_in1.dist.params[:m]*msg_in1.dist.params[:m]')
end
ruleSPMultiplicationOutVGP(msg_out::Void, msg_in1::Message{Gaussian, Univariate}, msg_a::Message{PointMass, Multivariate}) = ruleSPMultiplicationOutVPG(nothing, msg_a, msg_in1)

function ruleSPMultiplicationOutVGP(msg_out::Void, 
                                    msg_in1::Message{Gaussian, Multivariate}, 
                                    msg_a::Message{PointMass, Univariate})

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(Multivariate, Gaussian, m=msg_in1.dist.params[:m]*msg_a.dist.params[:m], v=msg_in1.dist.params[:v]*msg_a.dist.params[:m]^2)
end
ruleSPMultiplicationOutVPG(msg_out::Void, msg_in1::Message{PointMass, Univariate}, msg_a::Message{Gaussian, Multivariate}) = ruleSPMultiplicationOutVGP(nothing, msg_a, msg_in1)

# Multivariate inproduct # TODO: move to inproduct node
# function ruleSPMultiplicationIn1GVP(msg_out::Message{Gaussian, Multivariate},
#                                     msg_in1::Void,
#                                     msg_a::Message{PointMass, Multivariate})

#     ensureParameters!(msg_out.dist, (:xi, :w))
#     a = msg_a.dist.params[:m]
#     Message(Univariate, Gaussian, xi=a'*msg_out.dist.params[:xi], w=a'*msg_out.dist.params[:w]*a)
# end
# ruleSPMultiplicationAGPV(msg_out::Message{Gaussian, Multivariate}, msg_in1::Message{PointMass, Multivariate}, msg_a::Void) = ruleSPMultiplicationIn1GVP(msg_out, nothing, msg_in1)

# MatrixVariate*multivariate (in1 and a do NOT commute)
function ruleSPMultiplicationOutVGP(msg_out::Void,
                                    msg_in1::Message{Gaussian, Multivariate},
                                    msg_a::Message{PointMass, MatrixVariate})

    ensureParameters!(msg_in1.dist, (:m, :v))

    Message(Multivariate, Gaussian, m=msg_a.dist.params[:m]*msg_in1.dist.params[:m], v=msg_a.dist.params[:m]*msg_in1.dist.params[:v]*msg_a.dist.params[:m]')
end

function ruleSPMultiplicationIn1GVP(msg_out::Message{Gaussian, Multivariate},
                                    msg_in1::Void,
                                    msg_a::Message{PointMass, MatrixVariate})

    ensureParameters!(msg_out.dist, (:xi, :w))
    a = msg_a.dist.params[:m]
    w_q = a'*msg_out.dist.params[:w]*a
    w_q = w_q + tiny*diageye(size(w_q)[1]) # Ensure precision is always invertible
    Message(Multivariate, Gaussian, xi=a'*msg_out.dist.params[:xi], w=w_q)
end

ruleSPMultiplicationOutVPP(msg_out::Void, msg_in1::Message{PointMass, Multivariate}, msg_a::Message{PointMass, MatrixVariate}) = Message(Multivariate, PointMass, m=msg_a.dist.params[:m]*msg_in1.dist.params[:m])
ruleSPMultiplicationIn1PVP(msg_out::Message{PointMass, Multivariate}, msg_in1::Void, msg_a::Message{PointMass, MatrixVariate}) = Message(Multivariate, PointMass, m=pinv(msg_a.dist.params[:m])*msg_out.dist.params[:m])