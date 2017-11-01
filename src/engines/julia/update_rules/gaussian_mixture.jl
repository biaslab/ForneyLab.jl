export
ruleVBGaussianMixtureM1, 
ruleVBGaussianMixtureW1, 
ruleVBGaussianMixtureM2, 
ruleVBGaussianMixtureW2, 
ruleVBGaussianMixtureZ, 
ruleVBGaussianMixtureOut

function GMBackwardMRule(q_w_k::ProbabilityDistribution{Univariate}, q_x::ProbabilityDistribution{Univariate}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(Univariate(Gaussian, m=unsafeMean(q_x), w=z_k_hat*unsafeMean(q_w_k)))
end

function GMBackwardWRule(q_m_k::ProbabilityDistribution{Univariate}, q_x::ProbabilityDistribution{Univariate}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(Univariate(Gamma,
        a = 1.0 + 0.5*z_k_hat,
        b = 0.5*z_k_hat*( (unsafeMean(q_x) - unsafeMean(q_m_k))^2 + unsafeCov(q_x) + unsafeCov(q_m_k) ) ))
end

function GMBackwardZRule(   q_m::Vector,
                            q_w::Vector,
                            q_x::ProbabilityDistribution{Univariate})
    rho = zeros(2)
    for k = 1:2
        rho[k] = clamp(exp(-averageEnergy(GaussianMeanPrecision, q_x, q_m[k], q_w[k])), tiny, huge)
    end
    
    Message(Univariate(Bernoulli, p=rho[1]/sum(rho)))
end

function GMForwardRule( q_m::Vector,
                        q_w::Vector,
                        z_hat::Vector{Float64})
    w  = 0.0
    xi = 0.0
    for k = 1:length(z_hat)
        w  += z_hat[k]*unsafeMean(q_w[k])
        xi += unsafeMean(q_w[k])*unsafeMean(q_m[k])*z_hat[k]
    end

    Message(Univariate(Gaussian, xi=xi, w=w))
end

function ruleVBGaussianMixtureM1(   dist_out::ProbabilityDistribution{Univariate},
                                    dist_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    dist_mean_1::Any,
                                    dist_prec_1::ProbabilityDistribution{Univariate},
                                    dist_mean_2::ProbabilityDistribution{Univariate},
                                    dist_prec_2::ProbabilityDistribution{Univariate})

    GMBackwardMRule(dist_prec_1, dist_out, unsafeMean(dist_switch))
end

function ruleVBGaussianMixtureW1(   dist_out::ProbabilityDistribution{Univariate},
                                    dist_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    dist_mean_1::ProbabilityDistribution{Univariate},
                                    dist_prec_1::Any,
                                    dist_mean_2::ProbabilityDistribution{Univariate},
                                    dist_prec_2::ProbabilityDistribution{Univariate})

    GMBackwardWRule(dist_mean_1, dist_out, unsafeMean(dist_switch))
end

function ruleVBGaussianMixtureM2(   dist_out::ProbabilityDistribution{Univariate},
                                    dist_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    dist_mean_1::ProbabilityDistribution{Univariate},
                                    dist_prec_1::ProbabilityDistribution{Univariate},
                                    dist_mean_2::Any,
                                    dist_prec_2::ProbabilityDistribution{Univariate})

    GMBackwardMRule(dist_prec_2, dist_out, 1.0 - unsafeMean(dist_switch))
end

function ruleVBGaussianMixtureW2(   dist_out::ProbabilityDistribution{Univariate},
                                    dist_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    dist_mean_1::ProbabilityDistribution{Univariate},
                                    dist_prec_1::ProbabilityDistribution{Univariate},
                                    dist_mean_2::ProbabilityDistribution{Univariate},
                                    dist_prec_2::Any)

    GMBackwardWRule(dist_mean_2, dist_out, 1.0 - unsafeMean(dist_switch))
end

function ruleVBGaussianMixtureZ(dist_out::ProbabilityDistribution{Univariate},
                                dist_switch::Any,
                                dist_mean_1::ProbabilityDistribution{Univariate},
                                dist_prec_1::ProbabilityDistribution{Univariate},
                                dist_mean_2::ProbabilityDistribution{Univariate},
                                dist_prec_2::ProbabilityDistribution{Univariate})

    GMBackwardZRule([dist_mean_1, dist_mean_2], [dist_prec_1, dist_prec_2], dist_out)
end

function ruleVBGaussianMixtureOut(  dist_out::Any,
                                    dist_switch::ProbabilityDistribution{Univariate, Bernoulli},
                                    dist_mean_1::ProbabilityDistribution{Univariate},
                                    dist_prec_1::ProbabilityDistribution{Univariate},
                                    dist_mean_2::ProbabilityDistribution{Univariate},
                                    dist_prec_2::ProbabilityDistribution{Univariate})

    GMForwardRule([dist_mean_1, dist_mean_2], [dist_prec_1, dist_prec_2], [unsafeMean(dist_switch), 1.0 - unsafeMean(dist_switch)])
end
    