export
ruleSPProbitOutNG,
ruleSPProbitIn1PN,
ruleEPProbitIn1BG,
ruleEPProbitIn1CG,
ruleEPProbitIn1PG

function ruleSPProbitOutNG(msg_out::Nothing,
                           msg_in1::Message{F, Univariate}) where F<:Gaussian

    d_in1 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)
    
    p = normcdf(d_in1.params[:m] / sqrt(1 + d_in1.params[:v]))
    isnan(p) && (p = 0.5)

    Message(Univariate, Bernoulli, p=p)
end

function ruleSPProbitIn1PN(msg_out::Message{PointMass, Univariate},
                           msg_in1::Nothing)

    p = msg_out.dist.params[:m]
    if p == 1.0
        log_pdf = normlogcdf
    elseif p == 0.0
        log_pdf = (z)->normlogcdf(-z)
    else
        log_pdf = (z)->log(p*normcdf(z) + (1.0-p)*normcdf(-z))
    end

    return Message(Univariate, Function, log_pdf=log_pdf)
end

function ruleEPProbitIn1BG(msg_out::Message{Bernoulli, Univariate}, 
                           msg_in1::Message{F, Univariate}) where F<:Gaussian

    # Calculate approximate (Gaussian) message towards i[:in1]
    # The approximate message is an 'expectation' under the context (cavity distribution) encoded by incoming message msg_cavity.
    # Propagating the resulting approximate msg through the factor graph results in the expectation propagation (EP) algorithm.
    # Approximation procedure:
    #  1. Calculate exact (non-Gaussian) message towards i[:in1].
    #  2. Combine exact outbound msg on i[:in1] with exact inbound msg (cavity distribution) to find exact marginal.
    #  3. Approximate the exact (non-Gaussian) marginal with a Gaussian one using moment matching, under the constraint VAR[approximate marginal] > VAR[cavity].
    #  4. Calculate back the Gaussian outbound msg on i[:in1] that yields this approximate Gaussian marginal.
    # IMPORTANT NOTES:
    #  - This calculation results in an implicit cycle in the factor graph since the outbound message depends on the inbound message (cavity dist.).

    # Shorthand notations
    p = msg_out.dist.params[:p]
    (μ, σ2) = unsafeMeanCov(msg_in1.dist) # Moments of cavity distribution
    (ξ, γ) = unsafeWeightedMeanPrecision(msg_in1.dist)

    # Calculate first and second moment (mp_1, mp_2) of the 'true' marginal p(x) on edge connected to i[:in1]
    # p(x) = f(x) / Z
    # f(x) = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(x)*N(x|μ,σ2)
    #      = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(z)*(Φ(x)*N(x|μ,σ2)/Φ(z))
    #      = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(z)*g(x)
    # See paper for detailed derivation

    z = μ / sqrt(1 + σ2)
    N = exp(-0.5*z^2)./sqrt(2*pi) # 𝓝(z)

    # Moments of g(x)
    mg_1 = normcdf(z)*μ + σ2*N / sqrt(1+σ2)  # First moment of g
    mg_2 = 2*μ*mg_1 + normcdf(z)*(σ2 - μ^2) - σ2^2*z*N / (1+σ2)  # Second moment of g

    # Moments of f(x) (exact marginal)
    Z = 1 - p + (2*p-1)*normcdf(z)
    mp_1 = ((1-p)*μ + (2*p-1)*mg_1) / Z
    mp_2 = ((1-p)*(μ^2+σ2) + (2*p-1)*mg_2) / Z

    # Calculate Gaussian marginal with identical first and second moments (moment matching approximation)
    marginal_v = mp_2 - mp_1^2
    marginal_v = clamp(marginal_v, tiny, σ2-tiny) # ensure variance of marginal is not larger than variance of cavity distribution
    marginal_w = 1.0 / marginal_v
    marginal_xi = marginal_w * mp_1

    # Calculate the approximate message towards i[:in1]
    outbound_dist_w = marginal_w - γ
    outbound_dist_xi = marginal_xi - ξ

    return Message(Univariate, GaussianWeightedMeanPrecision, xi=outbound_dist_xi, w=outbound_dist_w)
end

function ruleEPProbitIn1PG(msg_out::Message{PointMass, Univariate}, msg_in1::Message{F, Univariate}) where F<:Gaussian
    p = msg_out.dist.params[:m]
    isnan(p) && (p = 0.5)
    (0.0 <= p <= 1.0) || error("Binary input $p must be between 0 and 1")

    return ruleEPProbitIn1BG(Message(Univariate, Bernoulli, p=p), msg_in1)
end

function ruleEPProbitIn1CG(msg_out::Message{Categorical, Univariate}, msg_in1::Message{F, Univariate}) where F<:Gaussian
    (length(msg_out.dist.params[:p]) == 2) || error("Probit node only supports categorical messages with 2 categories")
    p = msg_out.dist.params[:p][1]

    return ruleEPProbitIn1BG(Message(Univariate, Bernoulli, p=p), msg_in1)
end
