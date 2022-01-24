export
ruleSPGaussianPrecisionOutNPP,
ruleSPGaussianPrecisionMPNP,
ruleSPGaussianPrecisionOutNGP, 
ruleSPGaussianPrecisionMGNP, 
ruleVBGaussianPrecisionM, 
ruleVBGaussianPrecisionW, 
ruleVBGaussianPrecisionOut,
ruleSVBGaussianPrecisionOutVGD,
ruleSVBGaussianPrecisionW,
ruleSVBGaussianPrecisionMGVD,
ruleMGaussian{Precision}GGD,
ruleMGaussian{Precision}GGN

ruleSPGaussianPrecisionOutNPP(  msg_out::Nothing,
                                    msg_mean::Message{PointMass, V},
                                    msg_prec::Message{PointMass}) where V<:VariateType =
    Message(V, Gaussian{Precision}, m=deepcopy(msg_mean.dist.params[:m]), w=deepcopy(msg_prec.dist.params[:m]))

ruleSPGaussianPrecisionMPNP(msg_out::Message{PointMass}, msg_mean::Nothing, msg_prec::Message{PointMass}) = 
    ruleSPGaussianPrecisionOutNPP(msg_mean, msg_out, msg_prec)

function ruleSPGaussianPrecisionOutNGP( msg_out::Nothing,
                                            msg_mean::Message{F, V},
                                            msg_prec::Message{PointMass}) where {F<:Gaussian, V<:VariateType}

    d_mean = convert(Distribution{V, Gaussian{Moments}}, msg_mean.dist)

    Message(V, Gaussian{Moments}, m=d_mean.params[:m], v=d_mean.params[:v] + cholinv(msg_prec.dist.params[:m]))
end

ruleSPGaussianPrecisionMGNP(msg_out::Message{F}, msg_mean::Nothing, msg_prec::Message{PointMass}) where F<:Gaussian = 
    ruleSPGaussianPrecisionOutNGP(msg_mean, msg_out, msg_prec)

ruleVBGaussianPrecisionM(   dist_out::Distribution{V},
                                dist_mean::Any,
                                dist_prec::Distribution) where V<:VariateType =
    Message(V, Gaussian{Precision}, m=unsafeMean(dist_out), w=unsafeMean(dist_prec))

function ruleVBGaussianPrecisionW(  dist_out::Distribution{Univariate},
                                        dist_mean::Distribution{Univariate},
                                        dist_prec::Any)

    (m_mean, v_mean) = unsafeMeanCov(dist_mean)
    (m_out, v_out) = unsafeMeanCov(dist_out)

    Message(Univariate, Gamma, a=1.5, b=0.5*(v_mean + v_out + (m_mean - m_out)^2))
end

function ruleVBGaussianPrecisionW(  dist_out::Distribution{Multivariate},
                                        dist_mean::Distribution{Multivariate},
                                        dist_prec::Any)

    (m_mean, v_mean) = unsafeMeanCov(dist_mean)
    (m_out, v_out) = unsafeMeanCov(dist_out)

    Message(MatrixVariate, Wishart, v=cholinv( v_mean + v_out + (m_mean - m_out)*(m_mean - m_out)' ), nu=dims(dist_out)[1] + 2.0) 
end

ruleVBGaussianPrecisionOut( dist_out::Any,
                                dist_mean::Distribution{V},
                                dist_prec::Distribution) where V<:VariateType =
    Message(V, Gaussian{Precision}, m=unsafeMean(dist_mean), w=unsafeMean(dist_prec))

ruleSVBGaussianPrecisionOutVGD(dist_out::Any,
                                   msg_mean::Message{<:Gaussian, V},
                                   dist_prec::Distribution) where V<:VariateType = 
    Message(V, Gaussian{Moments}, m=unsafeMean(msg_mean.dist), v=unsafeCov(msg_mean.dist) + cholinv(unsafeMean(dist_prec)))

function ruleSVBGaussianPrecisionW(
    dist_out_mean::Distribution{Multivariate, F},
    dist_prec::Any) where F<:Gaussian

    joint_d = dims(dist_out_mean)[1]
    d_out_mean = convert(Distribution{Multivariate, Gaussian{Moments}}, dist_out_mean)
    (m, V) = unsafeMeanCov(d_out_mean)
    if joint_d == 2
        return Message(Univariate, Gamma, a=1.5, b=0.5*(V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2))
    else
        d = Int64(joint_d/2)
        return Message(MatrixVariate, Wishart, v=cholinv( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end, 1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ), nu=d + 2.0) 
    end
end

function ruleSVBGaussianPrecisionMGVD(  msg_out::Message{F, V},
                                            dist_mean::Any,
                                            dist_prec::Distribution) where {F<:Gaussian, V<:VariateType}

    d_out = convert(Distribution{V, Gaussian{Moments}}, msg_out.dist)

    Message(V, Gaussian{Moments}, m=d_out.params[:m], v=d_out.params[:v] + cholinv(unsafeMean(dist_prec)))
end

function ruleMGaussian{Precision}GGD(
    msg_out::Message{<:Gaussian, V},
    msg_mean::Message{<:Gaussian, V},
    dist_prec::Distribution) where V<:VariateType

    d_mean = convert(Distribution{V, Gaussian{Canonical}}, msg_mean.dist)
    d_out = convert(Distribution{V, Gaussian{Canonical}}, msg_out.dist)
    
    xi_y = d_out.params[:xi]
    W_y = d_out.params[:w]
    xi_m = d_mean.params[:xi]
    W_m = d_mean.params[:w]
    W_bar = unsafeMean(dist_prec)

    return Distribution(Multivariate, Gaussian{Canonical}, xi=[xi_y; xi_m], w=[W_y+W_bar -W_bar; -W_bar W_m+W_bar])
end

function ruleMGaussian{Precision}GGN(
    msg_out::Message{<:Gaussian, V},
    msg_mean::Message{<:Gaussian, V},
    msg_prec::Message{PointMass}) where V<:VariateType

    d_mean = convert(Distribution{V, Gaussian{Canonical}}, msg_mean.dist)
    d_out = convert(Distribution{V, Gaussian{Canonical}}, msg_out.dist)
    
    xi_y = d_out.params[:xi]
    W_y = d_out.params[:w]
    xi_m = d_mean.params[:xi]
    W_m = d_mean.params[:w]
    W_bar = msg_prec.dist.params[:m]

    return Distribution(Multivariate, Gaussian{Canonical}, xi=[xi_y; xi_m], w=[W_y+W_bar -W_bar; -W_bar W_m+W_bar])
end