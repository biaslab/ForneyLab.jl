############################################
# GaussianDistribution
############################################
# Description:
#   Encodes a univariate or multivariate Gaussian PDF.
#   Define (mean (m) or weighted mean (xi))
#   and (covariance (V) or precision (W)).
#   These result in the same PDF:
#    dist = GaussianDistribution(m=1.0, V=2.0)
#    dist = GaussianDistribution(m=1.0, W=0.5)
#    dist = GaussianDistribution(xi=0.5, W=0.5)
#    dist = GaussianDistribution(xi=0.5, V=2.0)
#   In case of a multivariate distribution,
#   m&xi should be vectors and V&W should be matrices.
############################################

export
    GaussianDistribution,
    ensureMVParametrization!,
    ensureMWParametrization!,
    ensureXiVParametrization!,
    ensureXiWParametrization!,
    isWellDefined,
    isConsistent

type GaussianDistribution <: ProbabilityDistribution
    m::Union(Vector{Float64}, Nothing)    # Mean vector
    V::Union(Matrix{Float64}, Nothing)       # Covariance matrix
    W::Union(Matrix{Float64}, Nothing)       # Weight matrix
    xi::Union(Vector{Float64}, Nothing)   # Weighted mean vector: xi=W*m
end
function GaussianDistribution(; m::Union(Float64,Vector{Float64},Nothing)=nothing,
                                V::Union(Float64,Matrix{Float64},Nothing)=nothing,
                                W::Union(Float64,Matrix{Float64},Nothing)=nothing,
                                xi::Union(Float64,Vector{Float64},Nothing)=nothing)
    m = (typeof(m)==Float64) ? [m] : deepcopy(m)
    V = (typeof(V)==Float64) ? fill!(Array(Float64,1,1),V) : deepcopy(V)
    W = (typeof(W)==Float64) ? fill!(Array(Float64,1,1),W) : deepcopy(W)
    xi = (typeof(xi)==Float64) ? [xi] : deepcopy(xi)

    self = GaussianDistribution(m, V, W, xi)

    # Check parameterizations
    isWellDefined(self) || error("Cannot create GaussianDistribution, parameterization is underdetermined.")

    return self
end
GaussianDistribution() = GaussianDistribution(m=0.0, V=1.0)
abstract BiVariateGaussianDistribution # Only used for uninformative function

uninformative(::Type{GaussianDistribution}) = GaussianDistribution(m=0.0, V=huge())
uninformative(::Type{BiVariateGaussianDistribution}) = GaussianDistribution(m=[0.0, 0.0], V=[huge() 0.0; 0.0 huge()])
uninformative(::Type{Float64}) = 1.0 # Float can be seen as a Gaussian with zero variance (delta peak at value)

function show(io::IO, dist::GaussianDistribution)
    println(io, "GaussianDistribution")
    print(io, "m  = ")
    show(io, dist.m)
    print(io, "\nV  = ")
    show(io, dist.V)
    print(io, "\nW  = ")
    show(io, dist.W)
    print(io, "\nxi = ")
    show(io, dist.xi)
    print(io, "\n")
end

Base.mean(dist::GaussianDistribution) = ensureMDefined!(dist).m
Base.var(dist::GaussianDistribution) = diag(ensureVDefined!(dist).V, 0)

# Methods to check and convert different parametrizations
function isWellDefined(dist::GaussianDistribution)
    # Check if dist is not underdetermined
    if ((is(dist.m, nothing) && is(dist.xi, nothing)) ||
        (is(dist.V, nothing) && is(dist.W, nothing)))
        return false
    end
    dimensions=0
    for field in [:m, :xi, :V, :W]
        if !is(getfield(dist, field), nothing)
            if dimensions>0
                if maximum(size(getfield(dist, field)))!=dimensions
                    return false
                end
            else
                dimensions = maximum(size(getfield(dist, field)))
            end
        end
    end
    return true
end
function isConsistent(dist::GaussianDistribution)
    # Check if dist is consistent in case it is overdetermined
    if !is(dist.V, nothing) && !is(dist.W, nothing)
        V_W_consistent = false
        try
           V_W_consistent = isApproxEqual(inv(dist.V), dist.W)
        catch
            try
                V_W_consistent = isApproxEqual(inv(dist.W), dist.V)
            catch
                error("Cannot check consistency of GaussianDistribution because both V and W are non-invertible.")
            end
        end
        if !V_W_consistent
            return false # V and W are not consistent
        end
    end
    if !is(dist.m, nothing) && !is(dist.xi, nothing)
        if !is(dist.V, nothing)
            if isApproxEqual(dist.V * dist.xi, dist.m) == false
                return false # m and xi are not consistent
            end
        else
            if isApproxEqual(dist.W * dist.m, dist.xi) == false
                return false # m and xi are not consistent
            end
        end
    end
    return true # all validations passed
end
function ensureMDefined!(dist::GaussianDistribution)
    # Ensure that dist.m is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    dist.m = is(dist.m, nothing) ? ensureVDefined!(dist).V * dist.xi : dist.m
    return dist
end
function ensureXiDefined!(dist::GaussianDistribution)
    # Ensure that dist.xi is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    dist.xi = is(dist.xi, nothing) ? ensureWDefined!(dist).W * dist.m : dist.xi
    return dist
end
function ensureVDefined!(dist::GaussianDistribution)
    # Ensure that dist.V is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    try
        dist.V = is(dist.V, nothing) ? inv(dist.W) : dist.V
    catch
        error("Cannot calculate V of GaussianDistribution because W is not invertible.")
    end
    return dist
end
function ensureWDefined!(dist::GaussianDistribution)
    # Ensure that dist.W is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    try
        dist.W = is(dist.W, nothing) ? inv(dist.V) : dist.W
    catch
        error("Cannot calculate W of GaussianDistribution because V is not invertible.")
    end
    return dist
end
ensureMVParametrization!(dist::GaussianDistribution) = ensureVDefined!(ensureMDefined!(dist))
ensureMWParametrization!(dist::GaussianDistribution) = ensureWDefined!(ensureMDefined!(dist))
ensureXiVParametrization!(dist::GaussianDistribution) = ensureVDefined!(ensureXiDefined!(dist))
ensureXiWParametrization!(dist::GaussianDistribution) = ensureWDefined!(ensureXiDefined!(dist))

function ==(x::GaussianDistribution, y::GaussianDistribution)
    if is(x, y) return true end
    eps = tiny()
    # Check m or xi
    if x.m!=nothing && y.m!=nothing
        (length(x.m)==length(x.m)) || return false
        (maximum(abs(x.m-y.m)) < eps) || return false
    elseif x.xi!=nothing && y.xi!=nothing
        (length(x.xi)==length(x.xi)) || return false
        (maximum(abs(x.xi-y.xi)) < eps) || return false
    else
        ensureMDefined!(x); ensureMDefined!(y);
        (length(x.m)==length(x.m)) || return false
        (maximum(abs(x.m-y.m)) < eps) || return false
    end

    # Check V or W
    if x.V!=nothing && y.V!=nothing
        (length(x.V)==length(x.V)) || return false
        (maximum(abs(x.V-y.V)) < eps) || return false
    elseif x.W!=nothing && y.W!=nothing
        (length(x.W)==length(x.W)) || return false
        (maximum(abs(x.W-y.W)) < eps) || return false
    else
        ensureVDefined!(x); ensureVDefined!(y);
        (length(x.V)==length(x.V)) || return false
        (maximum(abs(x.V-y.V)) < eps) || return false
    end

    return true
end