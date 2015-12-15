############################################
# WishartDistribution
############################################
# Description:
#   Encodes a Wishart distribution over precision
#   matrix W(W | V, nu).
#   (Bishop, 2007; "Pattern recognition and machine learning").
#   
#   Define:
#       Scale matrix V (p x p positive definite)), and
#       degrees of freedom nu > p - 1.
#   
#   Example:
#       WishartDistribution(V = 1.0, V=[2.0, 0.0; 0.0, 2.0])
############################################

export
    WishartDistribution,
    isProper

type WishartDistribution <: MultivariateProbabilityDistribution
    V::Matrix{Float64}  # Scale matrix
    nu::Float64         # Degrees of freedom
end

WishartDistribution(; V::Matrix{Float64} = [1.0].', nu::Float64 = 1.0) = WishartDistribution(V, nu)
WishartDistribution() = WishartDistribution(V = [1.0].', nu = 1.0)

vague(::Type{WishartDistribution}; dim=1) = WishartDistribution(V = huge*eye(dim), nu = tiny) # Scale invariant (Jeffrey's) prior in each dimension

function format(dist::WishartDistribution)
    if isProper(dist)
        return "W(V=$(format(dist.V)), ν=$(format(dist.nu)))"
    else
        return "W(underdetermined)"
    end
end

show(io::IO, dist::WishartDistribution) = println(io, format(dist))

function Base.mean(dist::WishartDistribution)
    if isProper(dist)
        return dist.nu*dist.V
    else
        return fill!(similar(dist.V), NaN)
    end
end

function Base.cov(dist::WishartDistribution)
    d = size(dist.V, 1)
    M = fill!(similar(dist.V), NaN)
    if isProper(dist)
        for i = 1:d
            for j = 1:d
                M[i, j] = dist.nu*(dist.V[i, j]^2 + dist.V[i, i]*dist.V[j, j])
            end
        end
        return M
    else
        return M
    end
end

function Base.var(dist::WishartDistribution)
    if isProper(dist)
        return dist.nu*2.0*diag(dist.V).^2
    else
        return fill!(zeros(d), NaN)
    end
end

function isProper(dist::WishartDistribution)
    (size(dist.V, 1) == size(dist.V, 2)) || return false
    (dist.nu > size(dist.V, 1) - 1) || return false
    isRoundedPosDef(dist.V) || return false

    return true
end

function ==(x::WishartDistribution, y::WishartDistribution)
    is(x, y) && return true
    isApproxEqual(x.nu, y.nu) || return false
    (length(x.V)==length(y.V)) || return false
    isApproxEqual(x.V, y.V) || return false

    return true
end

# Convert GammaDistribution -> WishartDistribution
convert(::Type{WishartDistribution}, d::GammaDistribution) = WishartDistribution(V = [1.0/(2.0*d.b)].', nu = 2.0*d.a)

# Convert WishartDistribution -> GammaDistribution
function convert(::Type{GammaDistribution}, d::WishartDistribution)
    (length(d.V) == 1) || error("Can only convert WishartDistribution to GammaDistribution if it has dimensionality 1")
    GammaDistribution(a = d.nu/2.0, b = 1.0/(2.0*d.V[1, 1]))
end
