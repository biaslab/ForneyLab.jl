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
ruleVBGaussianMomentsOut

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

    Message(Univariate, Function, log_pdf=(x)-> -0.5*log(d_out.params[:v] + d_mean.params[:v] + x) - 1/(2*x)*(d_out.params[:m] - d_mean.params[:m])^2)
end

function ruleSPGaussianMomentsVPGN(msg_out::Message{PointMass, Univariate},
                                   msg_mean::Message{F, Univariate},
                                   msg_var::Nothing) where F<:Gaussian

    d_mean = convert(Distribution{Univariate, Gaussian{Moments}}, msg_mean.dist)

    Message(Univariate, Function, log_pdf=(x)-> -0.5*log(d_mean.params[:v] + x) - 1/(2*x)*(msg_out.dist.params[:m] - d_mean.params[:m])^2)
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