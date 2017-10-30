export
ruleSPMultiplicationOutGP,
ruleSPMultiplicationInGP

function ruleSPMultiplicationOutGP( msg_out::Void,
                                    msg_in::Message{Univariate{Gaussian}},
                                    msg_gain::Message{Univariate{PointMass}})

    ensureParameters!(msg_in.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_gain.dist.params[:m]*msg_in.dist.params[:m], v=msg_gain.dist.params[:m]^2 * msg_in.dist.params[:v]))
end

function ruleSPMultiplicationInGP(  msg_out::Message{Univariate{Gaussian}},
                                    msg_in::Void,
                                    msg_gain::Message{Univariate{PointMass}})

    ensureParameters!(msg_out.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_out.dist.params[:m]/msg_gain.dist.params[:m], v=msg_out.dist.params[:v]/msg_gain.dist.params[:m]^2))
end