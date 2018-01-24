export
ruleVBGaussianMixtureM, 
ruleVBGaussianMixtureW, 
ruleVBGaussianMixtureZBer, 
ruleVBGaussianMixtureZCat, 
ruleVBGaussianMixtureOut

function ruleVBGaussianMixtureM(dist_out::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution,
                                dist_factors::Vararg{Union{Void, ProbabilityDistribution{Univariate}}})
    # Univariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    k = findfirst(dist_means .== nothing) # Find factor
    z_bar = clamp.(unsafeMeanVector(dist_switch), tiny, 1.0 - tiny)

    return Message(Univariate, Gaussian, m=unsafeMean(dist_out), w=z_bar[k]*unsafeMean(dist_precs[k]))
end

function ruleVBGaussianMixtureM(dist_out::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution,
                                dist_factors::Vararg{Union{Void, ProbabilityDistribution{Multivariate}, ProbabilityDistribution{MatrixVariate}}})
    # Multivariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    k = findfirst(dist_means .== nothing) # Find factor
    z_bar = clamp.(unsafeMeanVector(dist_switch), tiny, 1.0 - tiny)

    return Message(Multivariate, Gaussian, m=unsafeMean(dist_out), w=z_bar[k]*unsafeMean(dist_precs[k]))
end

function ruleVBGaussianMixtureW(dist_out::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution,
                                dist_factors::Vararg{Union{Void, ProbabilityDistribution{Univariate}}})
    # Univariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    k = findfirst(dist_precs .== nothing) # Find factor
    z_bar = unsafeMeanVector(dist_switch)

    return Message(Univariate, Gamma,
        a = 1.0 + 0.5*z_bar[k],
        b = 0.5*z_bar[k]*( (unsafeMean(dist_out) - unsafeMean(dist_means[k]))^2 + unsafeCov(dist_out) + unsafeCov(dist_means[k]) ) )
end

function ruleVBGaussianMixtureW(dist_out::ProbabilityDistribution,
                                dist_switch::ProbabilityDistribution,
                                dist_factors::Vararg{Union{Void, ProbabilityDistribution{Multivariate}, ProbabilityDistribution{MatrixVariate}}})
    # Multivariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    k = findfirst(dist_precs .== nothing) # Find factor
    z_bar = unsafeMeanVector(dist_switch)
    d = dims(dist_means[1])

    return Message(MatrixVariate, Wishart,
        nu = 1.0 + z_bar[k] + d,
        v  = cholinv( z_bar[k]*( (unsafeMean(dist_out) - unsafeMean(dist_means[k]))*(unsafeMean(dist_out) - unsafeMean(dist_means[k]))' + unsafeCov(dist_out) + unsafeCov(dist_means[k]) ) ))
end

function softmax(v::Vector{Float64})
    r = v - maximum(v)
    clamp!(r, -100.0, 0.0)
    exp.(r)./sum(exp.(r))
end

function ruleVBGaussianMixtureZBer( dist_out::ProbabilityDistribution,
                                    dist_switch::Any,
                                    dist_m1::ProbabilityDistribution,
                                    dist_w1::ProbabilityDistribution,
                                    dist_m2::ProbabilityDistribution,
                                    dist_w2::ProbabilityDistribution)
    # Uni- and Multivariate update
    U = Vector{Float64}(2)
    U[1] = averageEnergy(GaussianMeanPrecision, dist_out, dist_m1, dist_w1)
    U[2] = averageEnergy(GaussianMeanPrecision, dist_out, dist_m2, dist_w2)

    return Message(Univariate, Bernoulli, p=softmax(-U)[1])
end

function ruleVBGaussianMixtureZCat( dist_out::ProbabilityDistribution,
                                    dist_switch::Any,
                                    dist_factors::Vararg{ProbabilityDistribution})
    # Uni- and Multivariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])

    n_factors = length(dist_means)
    U = Vector{Float64}(n_factors)
    for k = 1:n_factors
        U[k] = averageEnergy(GaussianMeanPrecision, dist_out, dist_means[k], dist_precs[k])
    end

    return Message(Univariate, Categorical, p=softmax(-U))
end

function ruleVBGaussianMixtureOut(  dist_out::Any,
                                    dist_switch::ProbabilityDistribution,
                                    dist_factors::Vararg{ProbabilityDistribution{Univariate}})
    # Univariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    z_bar = unsafeMeanVector(dist_switch)

    w  = 0.0
    xi = 0.0
    for k = 1:length(z_bar)
        w  += z_bar[k]*unsafeMean(dist_precs[k])
        xi += unsafeMean(dist_precs[k])*unsafeMean(dist_means[k])*z_bar[k]
    end

    return Message(Univariate, Gaussian, xi=xi, w=w)
end

function ruleVBGaussianMixtureOut(  dist_out::Any,
                                    dist_switch::ProbabilityDistribution,
                                    dist_factors::Vararg{Union{ProbabilityDistribution{Multivariate}, ProbabilityDistribution{MatrixVariate}}})
    # Multivariate update
    dist_means = collect(dist_factors[1:2:end])
    dist_precs = collect(dist_factors[2:2:end])
    z_bar = unsafeMeanVector(dist_switch)
    d = dims(dist_means[1])

    w = Diagonal(zeros(d))
    xi = zeros(d)
    for k = 1:length(z_bar)
        w  += z_bar[k]*unsafeMean(dist_precs[k])
        xi += unsafeMean(dist_precs[k])*unsafeMean(dist_means[k])*z_bar[k]
    end

    return Message(Multivariate, Gaussian, xi=xi, w=w)
end    