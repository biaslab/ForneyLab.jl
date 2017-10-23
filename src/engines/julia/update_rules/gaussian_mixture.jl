export ruleVBGaussianMixture1, ruleVBGaussianMixture2, ruleVBGaussianMixture3, ruleVBGaussianMixture4, ruleVBGaussianMixture5, ruleVBGaussianMixture6

function GMBackwardMRule(q_w_k::ProbabilityDistribution, q_x::ProbabilityDistribution, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(Gaussian, m=unsafeMean(q_x), w=z_k_hat*unsafeMean(q_w_k))
end

function GMBackwardWRule(q_m_k::ProbabilityDistribution, q_x::ProbabilityDistribution, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(Gamma,
            a = 1.0 + 0.5*z_k_hat,
            b = 0.5*z_k_hat*( (unsafeMean(q_x) - unsafeMean(q_m_k))^2 + unsafeCov(q_x) + unsafeCov(q_m_k) ))
end

function GMBackwardZRule(   q_m::Vector,
                            q_w::Vector,
                            q_x::ProbabilityDistribution)
    rho = zeros(2)
    for k = 1:2
        rho[k] = clamp(exp(-averageEnergy(GaussianMeanPrecision, q_m[k], q_w[k], q_x)), tiny, huge)
    end
    
    Message(Bernoulli, p=rho[1]/sum(rho))
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

    Message(Gaussian, xi=xi, w=w)
end

function ruleVBGaussianMixture1(dist_mean_1::Any,
                                dist_prec_1::ProbabilityDistribution,
                                dist_mean_2::ProbabilityDistribution,
                                dist_prec_2::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution{Bernoulli},
                                dist_out::ProbabilityDistribution)

    GMBackwardMRule(dist_prec_1, dist_out, unsafeMean(dist_switch))
end

function ruleVBGaussianMixture2(dist_mean_1::ProbabilityDistribution,
                                dist_prec_1::Any,
                                dist_mean_2::ProbabilityDistribution,
                                dist_prec_2::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution{Bernoulli},
                                dist_out::ProbabilityDistribution)

    GMBackwardWRule(dist_mean_1, dist_out, unsafeMean(dist_switch))
end

function ruleVBGaussianMixture3(dist_mean_1::ProbabilityDistribution,
                                dist_prec_1::ProbabilityDistribution,
                                dist_mean_2::Any,
                                dist_prec_2::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution{Bernoulli},
                                dist_out::ProbabilityDistribution)

    GMBackwardMRule(dist_prec_2, dist_out, 1.0 - unsafeMean(dist_switch))
end

function ruleVBGaussianMixture4(dist_mean_1::ProbabilityDistribution,
                                dist_prec_1::ProbabilityDistribution,
                                dist_mean_2::ProbabilityDistribution,
                                dist_prec_2::Any,
                                dist_switch::ProbabilityDistribution{Bernoulli},
                                dist_out::ProbabilityDistribution)

    GMBackwardWRule(dist_mean_2, dist_out, 1.0 - unsafeMean(dist_switch))
end

function ruleVBGaussianMixture5(dist_mean_1::ProbabilityDistribution,
                                dist_prec_1::ProbabilityDistribution,
                                dist_mean_2::ProbabilityDistribution,
                                dist_prec_2::ProbabilityDistribution,
                                dist_switch::Any,
                                dist_out::ProbabilityDistribution)

    GMBackwardZRule([dist_mean_1, dist_mean_2], [dist_prec_1, dist_prec_2], dist_out)
end

function ruleVBGaussianMixture6(dist_mean_1::ProbabilityDistribution,
                                dist_prec_1::ProbabilityDistribution,
                                dist_mean_2::ProbabilityDistribution,
                                dist_prec_2::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution{Bernoulli},
                                dist_out::Any)

    GMForwardRule([dist_mean_1, dist_mean_2], [dist_prec_1, dist_prec_2], [unsafeMean(dist_switch), 1.0 - unsafeMean(dist_switch)])
end
    