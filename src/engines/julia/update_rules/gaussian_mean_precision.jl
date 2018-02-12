export
ruleSPGaussianMeanPrecisionOutVPP,
ruleSPGaussianMeanPrecisionMPVP,
ruleSPGaussianMeanPrecisionOutVGP, 
ruleSPGaussianMeanPrecisionMGVP, 
ruleVBGaussianMeanPrecisionM, 
ruleVBGaussianMeanPrecisionW, 
ruleVBGaussianMeanPrecisionOut,
ruleSVBGaussianMeanPrecisionOutVGD,
ruleSVBGaussianMeanPrecisionW,
ruleSVBGaussianMeanPrecisionMGVD,
ruleMGaussianMeanPrecisionGGD

ruleSPGaussianMeanPrecisionOutVPP{V<:VariateType}(  msg_out::Void,
                                                    msg_mean::Message{PointMass, V},
                                                    msg_prec::Message{PointMass}) =
    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), w=deepcopy(msg_prec.dist.params[:m]))

ruleSPGaussianMeanPrecisionMPVP(msg_out::Message{PointMass}, msg_mean::Void, msg_prec::Message{PointMass}) = ruleSPGaussianMeanPrecisionOutVPP(msg_mean, msg_out, msg_prec)

function ruleSPGaussianMeanPrecisionOutVGP{V<:VariateType}( msg_out::Void,
                                                            msg_mean::Message{Gaussian, V},
                                                            msg_prec::Message{PointMass})

    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), v=msg_mean.dist.params[:v] + cholinv(msg_prec.dist.params[:m]))
end

ruleSPGaussianMeanPrecisionMGVP(msg_out::Message{Gaussian}, msg_mean::Void, msg_prec::Message{PointMass}) = ruleSPGaussianMeanPrecisionOutVGP(msg_mean, msg_out, msg_prec)

ruleVBGaussianMeanPrecisionM{V<:VariateType}(   dist_out::ProbabilityDistribution{V},
                                                dist_mean::Any,
                                                dist_prec::ProbabilityDistribution) =
    Message(V, Gaussian, m=unsafeMean(dist_out), w=unsafeMean(dist_prec))

ruleVBGaussianMeanPrecisionW(   dist_out::ProbabilityDistribution{Univariate},
                                dist_mean::ProbabilityDistribution{Univariate},
                                dist_prec::Any) =
    Message(Univariate, Gamma, a=1.5, b=0.5*(unsafeVar(dist_mean) + unsafeVar(dist_out) + (unsafeMean(dist_mean) - unsafeMean(dist_out))^2))

ruleVBGaussianMeanPrecisionW(   dist_out::ProbabilityDistribution{Multivariate},
                                dist_mean::ProbabilityDistribution{Multivariate},
                                dist_prec::Any) =
    Message(MatrixVariate, Wishart, v=cholinv( unsafeCov(dist_out) + unsafeCov(dist_mean) + (unsafeMean(dist_out) - unsafeMean(dist_mean))*(unsafeMean(dist_out) - unsafeMean(dist_mean))' ), nu=dims(dist_out) + 2.0) 

ruleVBGaussianMeanPrecisionOut{V<:VariateType}( dist_out::Any,
                                                dist_mean::ProbabilityDistribution{V},
                                                dist_prec::ProbabilityDistribution) =
    Message(V, Gaussian, m=unsafeMean(dist_mean), w=unsafeMean(dist_prec))

function ruleSVBGaussianMeanPrecisionOutVGD{V<:VariateType}(dist_out::Any,
                                                            msg_mean::Message{Gaussian, V},
                                                            dist_prec::ProbabilityDistribution)
    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(V, Gaussian, m=deepcopy(msg_mean.dist.params[:m]), v=msg_mean.dist.params[:v] + cholinv(unsafeMean(dist_prec)))
end

function ruleSVBGaussianMeanPrecisionW( dist_out_mean::ProbabilityDistribution{Multivariate, Gaussian},
                                        dist_prec::Any)
    ensureParameters!(dist_out_mean, (:m, :v))

    joint_dims = dims(dist_out_mean)
    if joint_dims == 2
        V = dist_out_mean.params[:v]
        m = dist_out_mean.params[:m]
        
        return Message(Univariate, Gamma, a=1.5, b=0.5*(V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2))
    else
        V = dist_out_mean.params[:v]
        m = dist_out_mean.params[:m]
        d = Int64(joint_dims/2)
        
        return Message(MatrixVariate, Wishart, v=cholinv( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end, 1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ), nu=d + 2.0) 
    end
end

function ruleSVBGaussianMeanPrecisionMGVD{V<:VariateType}(  msg_out::Message{Gaussian, V},
                                                            dist_mean::Any,
                                                            dist_prec::ProbabilityDistribution)
    ensureParameters!(msg_out.dist, (:m, :v))

    Message(V, Gaussian, m=deepcopy(msg_out.dist.params[:m]), v=msg_out.dist.params[:v] + cholinv(unsafeMean(dist_prec)))
end

function ruleMGaussianMeanPrecisionGGD{V<:VariateType}( msg_out::Message{Gaussian, V},
                                                        msg_mean::Message{Gaussian, V},
                                                        dist_prec::ProbabilityDistribution)
    ensureParameters!(msg_out.dist, (:m, :w))
    ensureParameters!(msg_mean.dist, (:m, :w))

    m_y = msg_out.dist.params[:m]
    W_y = msg_out.dist.params[:w]
    m_m = msg_mean.dist.params[:m]
    W_m = msg_mean.dist.params[:w]
    W_bar = unsafeMean(dist_prec)

    V_q = cholinv([W_y+W_bar -W_bar; -W_bar W_m+W_bar])

    return ProbabilityDistribution(Multivariate, Gaussian, m=V_q*[W_y*m_y; W_m*m_m], v=V_q)
end