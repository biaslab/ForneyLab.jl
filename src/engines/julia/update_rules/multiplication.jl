export
ruleSPMultiplicationOutVPG,
ruleSPMultiplicationOutVGP,
ruleSPMultiplicationOutVPP,
ruleSPMultiplicationIn1GVP,
ruleSPMultiplicationIn1PVP,
ruleSPMultiplicationIn2GPV,
ruleSPMultiplicationIn2PPV

# Univariate
function ruleSPMultiplicationOutVPG(msg_out::Void,
                                    msg_in1::Message{PointMass, Univariate},
                                    msg_in2::Message{Gaussian, Univariate})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_in1.dist.params[:m]*msg_in2.dist.params[:m], v=msg_in2.dist.params[:v]*msg_in1.dist.params[:m]^2)
end
ruleSPMultiplicationOutVGP(msg_out::Void, msg_in1::Message{Gaussian, Univariate}, msg_in2::Message{PointMass, Univariate}) = ruleSPMultiplicationOutVPG(nothing, msg_in2, msg_in1)

function ruleSPMultiplicationIn1GVP(msg_out::Message{Gaussian, Univariate},
                                    msg_in1::Void,
                                    msg_in2::Message{PointMass, Univariate})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate, Gaussian, m=msg_out.dist.params[:m]/msg_in2.dist.params[:m], v=msg_out.dist.params[:v]/msg_in2.dist.params[:m]^2)
end
ruleSPMultiplicationIn2GPV(msg_out::Message{Gaussian, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_in2::Void) = ruleSPMultiplicationIn1GVP(msg_out, nothing, msg_in1)

ruleSPMultiplicationOutVPP(msg_out::Void, msg_in1::Message{PointMass, Univariate}, msg_in2::Message{PointMass, Univariate}) = Message(Univariate, PointMass, m=msg_in1.dist.params[:m]*msg_in2.dist.params[:m])
ruleSPMultiplicationIn1PVP(msg_out::Message{PointMass, Univariate}, msg_in1::Void, msg_in2::Message{PointMass, Univariate}) = Message(Univariate, PointMass, m=msg_out.dist.params[:m]/msg_in2.dist.params[:m])
ruleSPMultiplicationIn2PPV(msg_out::Message{PointMass, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_in2::Void) = Message(Univariate, PointMass, m=msg_out.dist.params[:m]/msg_in1.dist.params[:m])

# Multivariate
function ruleSPMultiplicationOutVPG(msg_out::Void,
                                    msg_in1::Message{PointMass, MatrixVariate},
                                    msg_in2::Message{Gaussian, Multivariate})

    ensureParameters!(msg_in2.dist, (:m, :v))

    Message(Multivariate, Gaussian, m=msg_in1.dist.params[:m]*msg_in2.dist.params[:m], v=msg_in1.dist.params[:m]*msg_in2.dist.params[:v]*msg_in1.dist.params[:m]')
end
ruleSPMultiplicationOutVGP(msg_out::Void, msg_in1::Message{Gaussian, Multivariate}, msg_in2::Message{PointMass, MatrixVariate}) = ruleSPMultiplicationOutVPG(nothing, msg_in2, msg_in1)

function ruleSPMultiplicationIn1GVP(msg_out::Message{Gaussian, Multivariate},
                                    msg_in1::Void,
                                    msg_in2::Message{PointMass, MatrixVariate})

    ensureParameters!(msg_out.dist, (:m, :v))
    A_inv = pinv(msg_in2.dist.params[:m])

    Message(Multivariate, Gaussian, m=A_inv*msg_out.dist.params[:m], v=A_inv*msg_out.dist.params[:v]*A_inv')
end
ruleSPMultiplicationIn2GPV(msg_out::Message{Gaussian, Multivariate}, msg_in1::Message{PointMass, MatrixVariate}, msg_in2::Void) = ruleSPMultiplicationIn1GVP(msg_out, nothing, msg_in1)

ruleSPMultiplicationOutVPP(msg_out::Void, msg_in1::Message{PointMass, MatrixVariate}, msg_in2::Message{PointMass, Multivariate}) = Message(Multivariate, PointMass, m=msg_in1.dist.params[:m]*msg_in2.dist.params[:m])
ruleSPMultiplicationOutVPP(msg_out::Void, msg_in1::Message{PointMass, Multivariate}, msg_in2::Message{PointMass, MatrixVariate}) = Message(Multivariate, PointMass, m=msg_in2.dist.params[:m]*msg_in1.dist.params[:m])
ruleSPMultiplicationIn1PVP(msg_out::Message{PointMass, Multivariate}, msg_in1::Void, msg_in2::Message{PointMass, MatrixVariate}) = Message(Multivariate, PointMass, m=pinv(msg_in2.dist.params[:m])*msg_out.dist.params[:m])
ruleSPMultiplicationIn2PPV(msg_out::Message{PointMass, Multivariate}, msg_in1::Message{PointMass, MatrixVariate}, msg_in2::Void) = Message(Multivariate, PointMass, m=pinv(msg_in1.dist.params[:m])*msg_out.dist.params[:m])
