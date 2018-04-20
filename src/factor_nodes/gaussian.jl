export Gaussian

abstract type Gaussian <: SoftFactor end

function prod!{F1<:Gaussian, F2<:Gaussian}(
    x::ProbabilityDistribution{Univariate, F1},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}=ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=1.0))

    z.params[:xi] = unsafeWeightedMean(x) + unsafeWeightedMean(y)
    z.params[:w] = unsafePrecision(x) + unsafePrecision(y)

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{Univariate, Gaussian},
                            y::ProbabilityDistribution{Univariate, PointMass},
                            z::ProbabilityDistribution{Univariate, PointMass}=ProbabilityDistribution(Univariate, PointMass, m=0.0))

    z.params[:m] = y.params[:m]
    return z
end

function prod!( x::ProbabilityDistribution{Multivariate, Gaussian},
                y::ProbabilityDistribution{Multivariate, Gaussian},
                z::ProbabilityDistribution{Multivariate, Gaussian}=ProbabilityDistribution(Multivariate, Gaussian, xi=[NaN], w=[NaN].'))

    ensureParameters!(x, (:xi, :w))
    ensureParameters!(y, (:xi, :w))

    invalidate!(z, :m)
    invalidate!(z, :v)
    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, Gaussian},
                            y::ProbabilityDistribution{Multivariate, PointMass},
                            z::ProbabilityDistribution{Multivariate, PointMass}=ProbabilityDistribution(Multivariate, PointMass, m=[NaN]))

    z.params[:m] = deepcopy(y.params[:m])
    return z
end