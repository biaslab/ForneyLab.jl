export
ruleVBGaussianMixtureM1, 
ruleVBGaussianMixtureW1, 
ruleVBGaussianMixtureM2, 
ruleVBGaussianMixtureW2, 
ruleVBGaussianMixtureZ, 
ruleVBGaussianMixtureOut

function GMBackwardMRule{T<:VariateType, U<:VariateType}(q_w_k::ProbabilityDistribution{T}, q_x::ProbabilityDistribution{U}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(U, Gaussian, m=unsafeMean(q_x), w=z_k_hat*unsafeMean(q_w_k))
end

function GMBackwardWRule(q_m_k::ProbabilityDistribution{Univariate}, q_x::ProbabilityDistribution{Univariate}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(Univariate, Gamma,
        a = 1.0 + 0.5*z_k_hat,
        b = 0.5*z_k_hat*( (unsafeMean(q_x) - unsafeMean(q_m_k))^2 + unsafeCov(q_x) + unsafeCov(q_m_k) ) )
end

function GMBackwardWRule(q_m_k::ProbabilityDistribution{Multivariate}, q_x::ProbabilityDistribution{Multivariate}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    Message(MatrixVariate, Wishart,
        nu = 1.0 + z_k_hat + dims(q_m_k),
        v  = cholinv( z_k_hat*( (unsafeMean(q_x) - unsafeMean(q_m_k))*(unsafeMean(q_x) - unsafeMean(q_m_k))' + unsafeCov(q_x) + unsafeCov(q_m_k) ) ))
end

function GMBackwardZRule(   q_m::Vector,
                            q_w::Vector,
                            q_x::ProbabilityDistribution)

    n_factors = length(q_m)
    (length(q_w) == n_factors) || error("Length of mean and precision input vectors must match")

    rho = zeros(n_factors)
    for k = 1:n_factors
        rho[k] = clamp(exp(-averageEnergy(GaussianMeanPrecision, q_x, q_m[k], q_w[k])), tiny, huge)
    end

    if n_factors == 2
        return Message(Univariate, Bernoulli, p=rho[1]/sum(rho))
    elseif n_factors > 2
        return Message(Univariate, Categorical, p=rho./sum(rho))
    else
        error("Number of factors must be 2 or greater")
    end
end

function GMForwardRule{T<:ProbabilityDistribution{Univariate}, U<:ProbabilityDistribution{Univariate}}( 
        q_m::Vector{T},
        q_w::Vector{U},
        z_hat::Vector{Float64})
    
    w  = 0.0
    xi = 0.0
    for k = 1:length(z_hat)
        w  += z_hat[k]*unsafeMean(q_w[k])
        xi += unsafeMean(q_w[k])*unsafeMean(q_m[k])*z_hat[k]
    end

    Message(Univariate, Gaussian, xi=xi, w=w)
end

function GMForwardRule{T<:ProbabilityDistribution{Multivariate}, U<:ProbabilityDistribution{MatrixVariate}}(
        q_m::Vector{T},
        q_w::Vector{U},
        z_hat::Vector{Float64})

    d = dims(q_m[1])
    w = Diagonal(zeros(d))
    xi = zeros(d)
    for k = 1:length(z_hat)
        w  += z_hat[k]*unsafeMean(q_w[k])
        xi += unsafeMean(q_w[k])*unsafeMean(q_m[k])*z_hat[k]
    end

    Message(Multivariate, Gaussian, xi=xi, w=w)
end

function ruleVBGaussianMixtureM1(   dist_out::ProbabilityDistribution,
                                    dist_switch::ProbabilityDistribution{Univariate},
                                    dist_factors::Vararg)

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    fact = findfirst(dist_means .== nothing)
    GMBackwardMRule(dist_precs[fact], dist_out, unsafeMeanVector(dist_switch)[fact])
end

function ruleVBGaussianMixtureW1(   dist_out::ProbabilityDistribution,
                                    dist_switch::ProbabilityDistribution{Univariate},
                                    dist_factors::Vararg)

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    fact = findfirst(dist_precs .== nothing)
    GMBackwardWRule(dist_means[fact], dist_out, unsafeMeanVector(dist_switch)[fact])
end

function ruleVBGaussianMixtureM2(   dist_out::ProbabilityDistribution,
                                    dist_switch::ProbabilityDistribution{Univariate},
                                    dist_factors::Vararg)

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    fact = findfirst(dist_means .== nothing)
    GMBackwardMRule(dist_precs[fact], dist_out, unsafeMeanVector(dist_switch)[fact])
end

function ruleVBGaussianMixtureW2(   dist_out::ProbabilityDistribution,
                                    dist_switch::ProbabilityDistribution{Univariate},
                                    dist_factors::Vararg)

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    fact = findfirst(dist_precs .== nothing)
    GMBackwardWRule(dist_means[fact], dist_out, unsafeMeanVector(dist_switch)[fact])
end

function ruleVBGaussianMixtureZ(dist_out::ProbabilityDistribution,
                                dist_switch::Any,
                                dist_factors::Vararg)

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    GMBackwardZRule(dist_means, dist_precs, dist_out)
end

function ruleVBGaussianMixtureOut(  dist_out::Any,
                                    dist_switch::ProbabilityDistribution{Univariate},
                                    dist_factors::Vararg)

    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    GMForwardRule(dist_means, dist_precs, unsafeMeanVector(dist_switch))
end
    