# Cummulative Gaussian (CDF of standard normal distribution)
Φ(x::Union{Float64, Vector{Float64}}) = 0.5*erfc(-x./sqrt(2.))

export
ruleSPSigmoidBinNG,
ruleEPSigmoidRealGB,
ruleEPSigmoidRealGC,
ruleEPSigmoidRealGP

function ruleSPSigmoidBinNG(msg_bin::Nothing,
                            msg_real::Message{F, Univariate}) where F<:Gaussian

    d_real = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_real.dist)
    
    p = Φ(d_real.params[:m] / sqrt(1 + d_real.params[:v]))
    isnan(p) && (p = 0.5)

    Message(Univariate, Bernoulli, p=p)
end

function ruleEPSigmoidRealGB(   msg_bin::Message{Bernoulli, Univariate}, 
                                msg_real::Message{F, Univariate}) where F<:Gaussian

    # Calculate approximate (Gaussian) message towards i[:real]
    # The approximate message is an 'expectation' under the context (cavity distribution) encoded by incoming message msg_cavity.
    # Propagating the resulting approximate msg through the factor graph results in the expectation propagation (EP) algorithm.
    # Approximation procedure:
    #  1. Calculate exact (non-Gaussian) message towards i[:real].
    #  2. Combine exact outbound msg on i[:real] with exact inbound msg (cavity distribution) to find exact marginal.
    #  3. Approximate the exact (non-Gaussian) marginal with a Gaussian one using moment matching, under the constraint VAR[approximate marginal] > VAR[cavity].
    #  4. Calculate back the Gaussian outbound msg on i[:real] that yields this approximate Gaussian marginal.
    # IMPORTANT NOTES:
    #  - This calculation results in an implicit cycle in the factor graph since the outbound message depends on the inbound message (cavity dist.).

    # Shorthand notations
    p = msg_bin.dist.params[:p]
    (μ, σ2) = unsafeMeanCov(msg_real.dist) # Moments of cavity distribution
    (ξ, γ) = unsafeWeightedMeanPrecision(msg_real.dist)

    # Calculate first and second moment (mp_1, mp_2) of the 'true' marginal p(x) on edge connected to i[:real]
    # p(x) = f(x) / Z
    # f(x) = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(x)*N(x|μ,σ2)
    #      = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(z)*(Φ(x)*N(x|μ,σ2)/Φ(z))
    #      = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(z)*g(x)
    # See paper for detailed derivation

    z = μ / sqrt(1 + σ2)
    N = exp(-0.5*z^2)./sqrt(2*pi) # 𝓝(z)

    # Moments of g(x)
    mg_1 = Φ(z)*μ + σ2*N / sqrt(1+σ2)  # First moment of g
    mg_2 = 2*μ*mg_1 + Φ(z)*(σ2 - μ^2) - σ2^2*z*N / (1+σ2)  # Second moment of g

    # Moments of f(x) (exact marginal)
    Z = 1 - p + (2*p-1)*Φ(z)
    mp_1 = ((1-p)*μ + (2*p-1)*mg_1) / Z
    mp_2 = ((1-p)*(μ^2+σ2) + (2*p-1)*mg_2) / Z

    # Calculate Gaussian marginal with identical first and second moments (moment matching approximation)
    marginal_v = mp_2 - mp_1^2
    marginal_v = clamp(marginal_v, tiny, σ2-tiny) # ensure variance of marginal is not larger than variance of cavity distribution
    marginal_w = 1.0 / marginal_v
    marginal_xi = marginal_w * mp_1

    # Calculate the approximate message towards i[:real]
    outbound_dist_w = marginal_w - γ
    outbound_dist_xi = marginal_xi - ξ

    return Message(Univariate, GaussianWeightedMeanPrecision, xi=outbound_dist_xi, w=outbound_dist_w)
end

function ruleEPSigmoidRealGP(msg_bin::Message{PointMass, Univariate}, msg_real::Message{F, Univariate}) where F<:Gaussian
    p = mapToBernoulliParameterRange(msg_bin.dist.params[:m])

    return ruleEPSigmoidRealGB(Message(Univariate, Bernoulli, p=p), msg_real)
end

function ruleEPSigmoidRealGC(msg_cat::Message{Categorical, Univariate}, msg_real::Message{F, Univariate}) where F<:Gaussian
    (length(msg_cat.dist.params[:p]) == 2) || error("Sigmoid node only supports categorical messages with 2 categories")
    p = msg_cat.dist.params[:p][1]

    return ruleEPSigmoidRealGB(Message(Univariate, Bernoulli, p=p), msg_real)
end
