export
FactorFunction,
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
struct ProbabilityDistribution{var_type<:VariateType, family<:FactorFunction}
    params::Dict
end

"""Sample multiple realizations from a probability distribution"""
sample(dist::ProbabilityDistribution, n_samples::Int64) = [sample(dist) for i in 1:n_samples] # TODO: individual samples can be optimized

"""Extract VariateType from dist"""
variateType(dist::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:FactorFunction} = V

show(io::IO, dist::ProbabilityDistribution) = println(io, format(dist))

matches(Ta::Type{Pa}, Tb::Type{Pb}) where {Pa<:ProbabilityDistribution, Pb<:ProbabilityDistribution} = (Pa<:Pb)
matches(::Type{Nothing}, ::Type{T}) where T<:ProbabilityDistribution = false

mean(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
mode(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMode(dist) : error("mode($(dist)) is undefined because the distribution is improper.")
var(dist::ProbabilityDistribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::ProbabilityDistribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")
logPdf(dist::ProbabilityDistribution) = isProper(dist) ? logPdf(dist) : error("logPdf($(dist)) is undefined.")

"""
`PointMass` is an abstract type used to describe point mass distributions.
It never occurs in a `FactorGraph`, but it is used as a probability distribution
type.
"""
abstract type PointMass <: DeltaFactor end

slug(::Type{PointMass}) = "δ"

format(dist::ProbabilityDistribution{V, PointMass}) where V<:VariateType = "$(slug(PointMass))(m=$(format(dist.params[:m])))"

dims(dist::ProbabilityDistribution{Univariate, PointMass}) = 1
dims(dist::ProbabilityDistribution{Multivariate, PointMass}) = length(dist.params[:m])
dims(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = size(dist.params[:m])

# PointMass distribution constructors
ProbabilityDistribution(::Type{Univariate}, ::Type{PointMass}; m::Number=1.0) = ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>m))
ProbabilityDistribution(::Type{PointMass}; m::Number=1.0) = ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>m))
ProbabilityDistribution(::Type{Multivariate}, ::Type{PointMass}; m::Vector=[1.0]) = ProbabilityDistribution{Multivariate, PointMass}(Dict(:m=>m))
ProbabilityDistribution(::Type{MatrixVariate}, ::Type{PointMass}; m::AbstractMatrix=transpose([1.0])) = ProbabilityDistribution{MatrixVariate, PointMass}(Dict(:m=>m))

unsafeMean(dist::ProbabilityDistribution{T, PointMass}) where T<:VariateType = deepcopy(dist.params[:m])

unsafeMode(dist::ProbabilityDistribution{T, PointMass}) where T<:VariateType = deepcopy(dist.params[:m])

unsafeMeanVector(dist::ProbabilityDistribution{Univariate, PointMass}) = [dist.params[:m]]
unsafeMeanVector(dist::ProbabilityDistribution{Multivariate, PointMass}) = deepcopy(dist.params[:m])

unsafeInverseMean(dist::ProbabilityDistribution{Univariate, PointMass}) = 1.0/dist.params[:m]
unsafeInverseMean(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = cholinv(dist.params[:m])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, PointMass}) = log(clamp(dist.params[:m], tiny, Inf))
unsafeLogMean(dist::ProbabilityDistribution{Multivariate, PointMass}) = log.(clamp.(dist.params[:m], tiny, Inf))
unsafeLogMean(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = log.(clamp.(dist.params[:m], tiny, Inf))

unsafeDetLogMean(dist::ProbabilityDistribution{Univariate, PointMass}) = log(dist.params[:m])
unsafeDetLogMean(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = log(det(dist.params[:m]))

unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, PointMass}) = log(1.0 - dist.params[:m])

unsafeVar(dist::ProbabilityDistribution{Univariate, PointMass}) = 0.0
unsafeVar(dist::ProbabilityDistribution{Multivariate, PointMass}) = zeros(dims(dist))

unsafeCov(dist::ProbabilityDistribution{Univariate, PointMass}) = 0.0
unsafeCov(dist::ProbabilityDistribution{Multivariate, PointMass}) = zeros(dims(dist), dims(dist))

unsafeMeanCov(dist::ProbabilityDistribution) = (unsafeMean(dist), unsafeCov(dist)) # Can be overloaded for efficiency

unsafeWeightedMeanPrecision(dist::ProbabilityDistribution) = (unsafeWeightedMean(dist), unsafePrecision(dist)) # Can be overloaded for efficiency

isProper(::ProbabilityDistribution{T, PointMass}) where T<:VariateType = true

# Probability distribution parametrized by function
slug(::Type{Function}) = "f"

format(dist::ProbabilityDistribution{V, Function}) where V<:VariateType = "$(dist.params)"

# Distribution constructors
ProbabilityDistribution(::Type{V}, ::Type{Function}; kwargs...) where V<:VariateType = ProbabilityDistribution{V, Function}(Dict{Symbol,Any}(kwargs))
ProbabilityDistribution(::Type{Function}; kwargs...) = ProbabilityDistribution{Univariate, Function}(Dict{Symbol,Any}(kwargs))

unsafeMode(dist::ProbabilityDistribution{T, Function}) where T<:VariateType = deepcopy(dist.params[:mode])

vague(::Type{Function}) = ProbabilityDistribution(Univariate, Function)

isProper(dist::ProbabilityDistribution{<:VariateType, Function}) = haskey(dist.params, :log_pdf)

logPdf(dist::ProbabilityDistribution{<:VariateType, Function}, x) = dist.params[:log_pdf](x)

# Convert VariateTypes
convert(::Type{ProbabilityDistribution{Multivariate, PointMass}}, dist::ProbabilityDistribution{Univariate, PointMass}) =
    ProbabilityDistribution(Multivariate, PointMass, m=[dist.params[:m]])
convert(::Type{ProbabilityDistribution{MatrixVariate, PointMass}}, dist::ProbabilityDistribution{Univariate, PointMass}) =
    ProbabilityDistribution(MatrixVariate, PointMass, m=mat(dist.params[:m]))
convert(::Type{ProbabilityDistribution{MatrixVariate, PointMass}}, dist::ProbabilityDistribution{Multivariate, PointMass}) =
    ProbabilityDistribution(MatrixVariate, PointMass, m=reshape(dist.params[:m], dims(dist), 1))

sample(dist::ProbabilityDistribution{T, PointMass}) where T<:VariateType = deepcopy(dist.params[:m])

"""
Compute conditional differential entropy: H(Y|X) = H(Y, X) - H(X)
"""
conditionalDifferentialEntropy(marg_joint::ProbabilityDistribution{Multivariate}, marg_condition::Vararg{ProbabilityDistribution}) = differentialEntropy(marg_joint) - sum([differentialEntropy(marg) for marg in marg_condition])

"""
This function ensures the argument expression is evaluated at runtime, allowing access to local variables
"""
pack(expr) = expr

function ==(t::ProbabilityDistribution{var_t, fam_t}, u::ProbabilityDistribution{var_u, fam_u}) where {fam_t, fam_u, var_t, var_u}
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
