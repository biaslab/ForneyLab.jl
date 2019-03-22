export GaussianWeightedMeanPrecision

abstract type GaussianWeightedMeanPrecision <: Gaussian end

slug(::Type{GaussianWeightedMeanPrecision}) = "𝒩"

format(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
# ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=Matrix{Float64}(I,1,1)) = ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>PDMat(w)))
function ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=Matrix{Float64}(I,1,1))

   # Check for type of input precision matrix
   if typeof(w) == Array{Float64,2}

       # If standard array, cast directly to full PDMat
       w = PDMat(w)

   elseif typeof(w) == Array{Array{Float64,2},2}

       # If 'mat' object (from helpers), cast (1,1) element to full PDMat
       w = PDMat(w[1])

   elseif typeof(w) == Diagonal{Float64,Array{Float64,1}}

       # If Diagional matrix, cast to PDiagMat
       w = PDiagMat(w)

   end
   return ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
end

dims(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = length(dist.params[:xi])

vague(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=PDMat(tiny*Matrix{Float64}(I,dims,dims)))

unsafeMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = inv(dist.params[:w])*dist.params[:xi]

unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = diag(inv(dist.params[:w]))

unsafeCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = inv(dist.params[:w])

function unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    v = inv(dist.params[:w])
    return (v*dist.params[:xi], v)
end

unsafeWeightedMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:xi]) # unsafe weighted mean

unsafePrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))

isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
isProper(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = isRoundedPosDef(dist.params[:w])

function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
    (t === u) && return true
    isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
    (t === u) && return true
    isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
end
