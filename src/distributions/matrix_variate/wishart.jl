export
    WishartDistribution,
    isProper

"""
Description:

    Encodes a Wishart distribution over precision matrix W(W | V, nu).
    (Bishop, 2007; 'Pattern recognition and machine learning').
   
Parameters:

    Scale matrix V (p x p positive definite)), and
    degrees of freedom nu > p - 1.
   
Construction:

    WishartDistribution(V=eye(3), nu=3.0)
    WishartDistribution(V=diageye(3), nu=3.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type WishartDistribution{dims} <: MatrixVariateProbabilityDistribution
    V::AbstractMatrix{Float64}  # Scale matrix
    nu::Float64                 # Degrees of freedom
end

WishartDistribution(; V::AbstractMatrix{Float64} = [1.0].', nu::Float64 = 1.0) = WishartDistribution{size(V, 1)}(deepcopy(V), nu)

WishartDistribution() = WishartDistribution(V = [1.0].', nu = 1.0)

function vague!{dims}(dist::WishartDistribution{dims})
    dist.V = huge*eye(dims)
    dist.nu = tiny
    return dist
end

vague{dims}(::Type{WishartDistribution{dims}}) = WishartDistribution(V=huge*diageye(dims), nu=tiny)

format(dist::WishartDistribution) = "W(V=$(format(dist.V)), ν=$(format(dist.nu)))"

show(io::IO, dist::WishartDistribution) = println(io, format(dist))

function Base.mean(dist::WishartDistribution)
    if isProper(dist)
        return dist.nu*dist.V
    else
        return fill!(similar(dist.V), NaN)
    end
end

function Base.var(dist::WishartDistribution)
    d = size(dist.V, 1)
    M = fill!(similar(Matrix(dist.V)), NaN)
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
function convert(::Type{GammaDistribution}, d::WishartDistribution{1})
    GammaDistribution(a = d.nu/2.0, b = 1.0/(2.0*d.V[1, 1]))
end
