export
ruleSPGaussianMomentsOutNPP,
ruleSPGaussianMomentsMPNP,
ruleSPGaussianMomentsOutNGP,
ruleSPGaussianMomentsMGNP,
ruleSPGaussianMomentsVGGN,
ruleSPGaussianMomentsVPGN,
ruleSPGaussianMomentsOutNSP,
ruleSPGaussianMomentsMSNP,
ruleSPGaussianMomentsOutNGS,
ruleSPGaussianMomentsMGNS,
ruleVBGaussianMomentsM,
ruleVBGaussianMomentsOut,
ruleSVBGaussianMomentsOutVGD,
ruleSVBGaussianMomentsMGVD,
ruleMGaussianMomentsGGD,
ruleMGaussianMomentsGGN

ruleSPGaussianMomentsOutNPP(msg_out::Nothing,
                            msg_mean::Message{PointMass, V},
                            msg_var::Message{PointMass}) where V<:VariateType =
    Message(V, Gaussian{Moments}, m=deepcopy(msg_mean.dist.params[:m]), v=deepcopy(msg_var.dist.params[:m]))

ruleSPGaussianMomentsMPNP(msg_out::Message{PointMass}, msg_mean::Nothing, msg_var::Message{PointMass}) =
    ruleSPGaussianMomentsOutNPP(msg_mean, msg_out, msg_var)

function ruleSPGaussianMomentsOutNGP(msg_out::Nothing,
                                     msg_mean::Message{F, V},
                                     msg_var::Message{PointMass}) where {F<:Gaussian, V<:VariateType}

    d_mean = convert(Distribution{V, Gaussian{Moments}}, msg_mean.dist)

    Message(V, Gaussian{Moments}, m=d_mean.params[:m], v=d_mean.params[:v] + msg_var.dist.params[:m])
end

ruleSPGaussianMomentsMGNP(msg_out::Message{F}, msg_mean::Nothing, msg_var::Message{PointMass}) where F<:Gaussian =
    ruleSPGaussianMomentsOutNGP(msg_mean, msg_out, msg_var)

function ruleSPGaussianMomentsVGGN(msg_out::Message{F1, Univariate},
                                   msg_mean::Message{F2, Univariate},
                                   msg_var::Nothing) where {F1<:Gaussian, F2<:Gaussian}

    d_out  = convert(Distribution{Univariate, Gaussian{Moments}}, msg_out.dist)
    d_mean = convert(Distribution{Univariate, Gaussian{Moments}}, msg_mean.dist)

    Message(Univariate, Function, log_pdf=(x)-> -0.5*log(d_out.params[:v] + d_mean.params[:v] + x) - 1/(2*(d_out.params[:v] + d_mean.params[:v] + x))*(d_out.params[:m] - d_mean.params[:m])^2)
end

function ruleSPGaussianMomentsVPGN(msg_out::Message{PointMass, Univariate},
                                   msg_mean::Message{F, Univariate},
                                   msg_var::Nothing) where F<:Gaussian

    d_mean = convert(Distribution{Univariate, Gaussian{Moments}}, msg_mean.dist)

    Message(Univariate, Function, log_pdf=(x)-> -0.5*log(d_mean.params[:v] + x) - 1/(2*(d_mean.params[:v] + x))*(msg_out.dist.params[:m] - d_mean.params[:m])^2)
end

# Particle update
function ruleSPGaussianMomentsOutNSP(msg_out::Nothing,
                                     msg_mean::Message{SampleList, V},
                                     msg_var::Message{PointMass}) where {V<:VariateType}

    samples = bootstrap(msg_mean.dist, msg_var.dist)
    weights = msg_mean.dist.params[:w]

    return Message(V, SampleList, s=samples, w=weights)
end

# Particle update
ruleSPGaussianMomentsMSNP(msg_out::Message{SampleList},
                          msg_mean::Nothing,
                          msg_var::Message{PointMass}) = ruleSPGaussianMomentsOutNSP(msg_mean, msg_out, msg_var)

# Particle update
function ruleSPGaussianMomentsOutNGS(msg_out::Nothing,
                                     msg_mean::Message{F, V},
                                     msg_var::Message{SampleList}) where {F<:Gaussian, V<:VariateType}

    samples = bootstrap(msg_mean.dist, msg_var.dist)
    weights = msg_var.dist.params[:w]

    Message(V, SampleList, s=samples, w=weights)
end

# Particle update
ruleSPGaussianMomentsMGNS(msg_out::Message{F},
                          msg_mean::Nothing,
                          msg_var::Message{SampleList}) where F<:Gaussian = ruleSPGaussianMomentsOutNGS(msg_mean, msg_out, msg_var)


ruleVBGaussianMomentsM(dist_out::Distribution{V},
                       dist_mean::Any,
                       dist_var::Distribution) where V<:VariateType =
    Message(V, Gaussian{Moments}, m=unsafeMean(dist_out), v=unsafeMean(dist_var))

ruleVBGaussianMomentsOut(dist_out::Any,
                         dist_mean::Distribution{V},
                         dist_var::Distribution) where V<:VariateType =
    Message(V, Gaussian{Moments}, m=unsafeMean(dist_mean), v=unsafeMean(dist_var))

ruleSVBGaussianMomentsOutVGD(dist_out::Any, # Only implemented for PointMass variance
                             msg_mean::Message{<:Gaussian, V},
                             dist_var::Distribution{<:VariateType, PointMass}) where V<:VariateType = 
    Message(V, Gaussian{Moments}, m=unsafeMean(msg_mean.dist), v=unsafeCov(msg_mean.dist) + unsafeMean(dist_var))

ruleSVBGaussianMomentsMGVD(msg_out::Message{F, V}, # Only implemented for PointMass variance
                           dist_mean::Any,
                           dist_var::Distribution{<:VariateType, PointMass}) where {F<:Gaussian, V<:VariateType} =
    Message(V, Gaussian{Moments}, m=unsafeMean(msg_out.dist), v=unsafeCov(msg_out.dist) + unsafeMean(dist_var))

function ruleMGaussianMomentsGGD( # Only implemented for PointMass variance
    msg_out::Message{<:Gaussian, V},
    msg_mean::Message{<:Gaussian, V},
    dist_var::Distribution{<:VariateType, PointMass}) where V<:VariateType

    d_mean = convert(Distribution{V, Gaussian{Canonical}}, msg_mean.dist)
    d_out = convert(Distribution{V, Gaussian{Canonical}}, msg_out.dist)
    
    xi_y = d_out.params[:xi]
    W_y = d_out.params[:w]
    xi_m = d_mean.params[:xi]
    W_m = d_mean.params[:w]
    W_bar = cholinv(unsafeMean(dist_var))

    return Distribution(Multivariate, Gaussian{Canonical}, xi=[xi_y; xi_m], w=[W_y+W_bar -W_bar; -W_bar W_m+W_bar])
end

function ruleMGaussianMomentsGGN(
    msg_out::Message{<:Gaussian, V},
    msg_mean::Message{<:Gaussian, V},
    msg_var::Message{PointMass}) where V<:VariateType

    d_mean = convert(Distribution{V, Gaussian{Canonical}}, msg_mean.dist)
    d_out = convert(Distribution{V, Gaussian{Canonical}}, msg_out.dist)
    
    xi_y = d_out.params[:xi]
    W_y = d_out.params[:w]
    xi_m = d_mean.params[:xi]
    W_m = d_mean.params[:w]
    W_bar = cholinv(msg_var.dist.params[:m])

    return Distribution(Multivariate, Gaussian{Canonical}, xi=[xi_y; xi_m], w=[W_y+W_bar -W_bar; -W_bar W_m+W_bar])
end