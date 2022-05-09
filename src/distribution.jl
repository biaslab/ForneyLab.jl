export
FactorFunction,
Distribution,
P,
ProbabilityDistribution,
Univariate,
Multivariate,
MatrixVariate,
PointMass,
mean,
mode,
var,
cov,
differentialEntropy,
conditionalDifferentialEntropy,
averageEnergy,
==,
vague,
sample,
dims,
logPdf

abstract type VariateType end
abstract type Univariate <: VariateType end
abstract type Multivariate <: VariateType end
abstract type MatrixVariate <: VariateType end

"""Types through which a probability distribution may be defined"""
const FactorFunction = Union{FactorNode, Function}

"""Encodes a probability distribution as a `FactorFunction` of type `family` with fixed interfaces"""
struct Distribution{var_type<:VariateType, family<:FactorFunction}
    params::Dict
end

"""Aliases for Distribution definition"""
const P = Distribution
const ProbabilityDistribution = Distribution # For backwards compatibility

"""Sample multiple realizations from a probability distribution"""
sample(dist::Distribution, n_samples::Int64) = [sample(dist) for i in 1:n_samples] # TODO: individual samples can be optimized

"""Extract VariateType from dist"""
variateType(::Distribution{V, <:FactorFunction}) where V<:VariateType = V

"""Extract VariateType from dims tuple"""
variateType(::Nothing) = Univariate # Default
variateType(::Tuple{}) = Univariate 
variateType(::Tuple{Int}) = Multivariate 
variateType(::Tuple{Int, Int}) = MatrixVariate

show(io::IO, dist::Distribution) = println(io, format(dist))

"""Parametric ineritance rule for Distribution, uses << symbol for ease of notation"""
# Parametric inheritance from nonparametric Distribution (where clauses are required for signature matching)
<<(::Type{Distribution}, ::Type{Distribution}) = true
<<(::Type{Distribution{Va}}, ::Type{Distribution}) where Va<:VariateType = true
<<(::Type{Distribution{Va, Fa}}, ::Type{Distribution}) where {Va<:VariateType, Fa<:FactorFunction} = true

# Parametric inheritance from Distribution{<:VariateType}
<<(::Type{Distribution{Va}}, ::Type{Distribution{Vb}}) where {Va<:VariateType, Vb<:VariateType} = (Va==Vb)
<<(::Type{Distribution{Va, Fa}}, ::Type{Distribution{Vb}}) where {Va<:VariateType, Vb<:VariateType, Fa<:FactorFunction} = (Va==Vb)

# Parametric inheritance from Distribution{<:VariateType, <:FactorFunction}
<<(::Type{Distribution{Va, Fa}}, ::Type{Distribution{Vb, Fb}}) where {Va<:VariateType, Vb<:VariateType, Fa<:FactorFunction, Fb<:FactorFunction} = (Va==Vb) && (Fa<:Fb)

mean(dist::Distribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
mode(dist::Distribution) = isProper(dist) ? unsafeMode(dist) : error("mode($(dist)) is undefined because the distribution is improper.")
var(dist::Distribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::Distribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")
logPdf(dist::Distribution) = isProper(dist) ? logPdf(dist) : error("logPdf($(dist)) is undefined.")

"""
`PointMass` is an abstract type used to describe point mass distributions.
It never occurs in a `FactorGraph`, but it is used as a probability distribution
type.
"""
abstract type PointMass <: DeltaFactor end

slug(::Type{PointMass}) = "δ"

format(dist::Distribution{V, PointMass}) where V<:VariateType = "$(slug(PointMass))(m=$(format(dist.params[:m])))"

dims(dist::Distribution{<:VariateType, PointMass}) = size(dist.params[:m])

# PointMass distribution constructors
Distribution(::Type{Univariate}, ::Type{PointMass}; m::Number=1.0) = Distribution{Univariate, PointMass}(Dict(:m=>m))
Distribution(::Type{PointMass}; m::Number=1.0) = Distribution{Univariate, PointMass}(Dict(:m=>m))
Distribution(::Type{Multivariate}, ::Type{PointMass}; m::Vector=[1.0]) = Distribution{Multivariate, PointMass}(Dict(:m=>m))
Distribution(::Type{MatrixVariate}, ::Type{PointMass}; m::AbstractMatrix=transpose([1.0])) = Distribution{MatrixVariate, PointMass}(Dict(:m=>m))

unsafeMean(dist::Distribution{T, PointMass}) where T<:VariateType = deepcopy(dist.params[:m])

unsafeMode(dist::Distribution{T, PointMass}) where T<:VariateType = deepcopy(dist.params[:m])

unsafeMeanVector(dist::Distribution{Univariate, PointMass}) = [dist.params[:m]]
unsafeMeanVector(dist::Distribution{Multivariate, PointMass}) = deepcopy(dist.params[:m])

unsafeInverseMean(dist::Distribution{Univariate, PointMass}) = 1.0/dist.params[:m]
unsafeInverseMean(dist::Distribution{MatrixVariate, PointMass}) = cholinv(dist.params[:m])

unsafeLogMean(dist::Distribution{Univariate, PointMass}) = log(clamp(dist.params[:m], tiny, Inf))
unsafeLogMean(dist::Distribution{Multivariate, PointMass}) = log.(clamp.(dist.params[:m], tiny, Inf))
unsafeLogMean(dist::Distribution{MatrixVariate, PointMass}) = log.(clamp.(dist.params[:m], tiny, Inf))

unsafeDetLogMean(dist::Distribution{Univariate, PointMass}) = log(dist.params[:m])
unsafeDetLogMean(dist::Distribution{MatrixVariate, PointMass}) = logdet(dist.params[:m])

unsafeMirroredLogMean(dist::Distribution{Univariate, PointMass}) = log(1.0 - dist.params[:m])

unsafeVar(dist::Distribution{Univariate, PointMass}) = 0.0
unsafeVar(dist::Distribution{Multivariate, PointMass}) = zeros(dims(dist)) # Vector

unsafeCov(dist::Distribution{Univariate, PointMass}) = 0.0
unsafeCov(dist::Distribution{Multivariate, PointMass}) = zeros(dims(dist)[1], dims(dist)[1]) # Matrix

unsafeMeanCov(dist::Distribution) = (unsafeMean(dist), unsafeCov(dist)) # Can be overloaded for efficiency

unsafeWeightedMeanPrecision(dist::Distribution) = (unsafeWeightedMean(dist), unsafePrecision(dist)) # Can be overloaded for efficiency

isProper(::Distribution{T, PointMass}) where T<:VariateType = true

# Probability distribution parametrized by function
slug(::Type{Function}) = "f"

format(dist::Distribution{V, Function}) where V<:VariateType = "$(dist.params)"

# Distribution constructors
Distribution(::Type{V}, ::Type{Function}; kwargs...) where V<:VariateType = Distribution{V, Function}(Dict{Symbol,Any}(kwargs))
Distribution(::Type{Function}; kwargs...) = Distribution{Univariate, Function}(Dict{Symbol,Any}(kwargs))

unsafeMode(dist::Distribution{T, Function}) where T<:VariateType = deepcopy(dist.params[:mode])

vague(::Type{Function}) = Distribution(Univariate, Function)

isProper(dist::Distribution{<:VariateType, Function}) = haskey(dist.params, :log_pdf)

logPdf(dist::Distribution{<:VariateType, Function}, x) = dist.params[:log_pdf](x)

# Convert VariateTypes
convert(::Type{Distribution{Multivariate, PointMass}}, dist::Distribution{Univariate, PointMass}) =
    Distribution(Multivariate, PointMass, m=[dist.params[:m]])
convert(::Type{Distribution{MatrixVariate, PointMass}}, dist::Distribution{Univariate, PointMass}) =
    Distribution(MatrixVariate, PointMass, m=mat(dist.params[:m]))
convert(::Type{Distribution{MatrixVariate, PointMass}}, dist::Distribution{Multivariate, PointMass}) =
    Distribution(MatrixVariate, PointMass, m=reshape(dist.params[:m], dims(dist)[1], 1))

sample(dist::Distribution{T, PointMass}) where T<:VariateType = deepcopy(dist.params[:m])

"""
Compute conditional differential entropy: H(Y|X) = H(Y, X) - H(X)
"""
conditionalDifferentialEntropy(marg_joint::Distribution{Multivariate}, marg_condition::Vararg{Distribution}) = differentialEntropy(marg_joint) - sum([differentialEntropy(marg) for marg in marg_condition])

"""
This function ensures the argument expression is evaluated at runtime, allowing access to local variables
"""
pack(expr) = expr

function ==(t::Distribution{var_t, fam_t}, u::Distribution{var_u, fam_u}) where {fam_t, fam_u, var_t, var_u}
    (fam_t == fam_u) || return false
    (var_t == var_u) || return false

    # Parameter values are checked for approximate equality
    t.params === u.params && return true
    if length(t.params) != length(u.params) return false end
    for pair in t.params
        if !in(pair, u.params, (a,b)->isapprox(a, b, atol=tiny))
            return false
        end
    end

    return true
end

function gaussianQuadrature(h::Function; m::Float64=0.0, v::Float64=1.0)
    abscissas = [-7.12581, -6.4095, -5.81223, -5.27555, -4.77716, -4.30555, -3.85376, -3.41717, -2.99249, -2.57725, -2.1695, -1.76765, -1.37038, -0.9765, -0.584979, -0.194841, 0.194841, 0.584979, 0.9765, 1.37038, 1.76765, 2.1695, 2.57725, 2.99249, 3.41717, 3.85376, 4.30555, 4.77716, 5.27555, 5.81223, 6.4095, 7.12581]
    weights = [7.31068e-23, 9.23174e-19, 1.19734e-15, 4.21501e-13, 5.93329e-11, 4.09883e-9, 1.57417e-7, 3.65059e-6, 5.41658e-5, 0.000536268, 0.00365489, 0.0175534, 0.0604581, 0.15127, 0.277458, 0.375238, 0.375238, 0.277458, 0.15127, 0.0604581, 0.0175534, 0.00365489, 0.000536268, 5.41658e-5, 3.65059e-6, 1.57417e-7, 4.09883e-9, 5.93329e-11, 4.21501e-13, 1.19734e-15, 9.23174e-19, 7.31068e-23]

    1/sqrt(pi) * sum([weights[i]*h(sqrt(2*v)*abscissas[i] + m) for i=1:32])
end
