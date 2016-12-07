export
    Wishart,
    isProper

"""
Description:

    Encodes a Wishart distribution over precision matrix W(W | V, nu).
    (Bishop, 2007; 'Pattern recognition and machine learning').

Parameters:

    Scale matrix V (p x p positive definite)), and
    degrees of freedom nu > p - 1.

Construction:

    Wishart(V=eye(3), nu=3.0)
    Wishart(V=diageye(3), nu=3.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type Wishart{dims} <: MatrixVariate{dims, dims}
    V::AbstractMatrix{Float64}  # Scale matrix
    nu::Float64                 # Degrees of freedom
end

Wishart(; V::AbstractMatrix{Float64} = [1.0].', nu::Float64 = 1.0) = Wishart{size(V, 1)}(deepcopy(V), nu)

function pdf(dist::Wishart, x::AbstractMatrix{Float64})
    error("TODO: pdf(::Wishart, x) is not yet implemented")
end

function vague!{dims}(dist::Wishart{dims})
    dist.V = huge*eye(dims)
    dist.nu = tiny
    return dist
end

vague{dims}(::Type{Wishart{dims}}) = Wishart(V=huge*diageye(dims), nu=tiny)

format(dist::Wishart) = "W(V=$(format(dist.V)), ν=$(format(dist.nu)))"

show(io::IO, dist::Wishart) = println(io, format(dist))

function prod!{dims}(x::Wishart{dims}, y::Wishart{dims}, z::Wishart{dims}=Wishart(V=eye(dims), nu=1.0))
    # Multiplication of 2 Wishart PDFs: p(z) = p(x) * p(y)
    # Derivation available in notebook
    z.V = x.V * cholinv(x.V + y.V) * y.V
    z.nu = x.nu + y.nu - dims - 1.0

    return z
end

@symmetrical function prod!{dims}(x::Wishart{dims}, y::MatrixDelta{Float64,dims,dims}, z::MatrixDelta{Float64,dims,dims}=deepcopy(y))
    # Multiplication of Wishart PDF with a delta
    (z.M == y.M) || (z.M[:] = y.M)

    return z
end

unsafeMean(dist::Wishart) = dist.nu*dist.V

function unsafeVar(dist::Wishart)
    d = size(dist.V, 1)
    M = fill!(similar(Matrix(dist.V)), NaN)
    for i = 1:d
        for j = 1:d
            M[i, j] = dist.nu*(dist.V[i, j]^2 + dist.V[i, i]*dist.V[j, j])
        end
    end
    return M
end

function unsafeDetLogMean{dims}(dist::Wishart{dims})
    sum([digamma(0.5*(dist.nu + 1 - i)) for i = 1:dims]) +
    dims*log(2) +
    log(det(dist.V))
end

function isProper(dist::Wishart)
    (size(dist.V, 1) == size(dist.V, 2)) || return false
    (dist.nu > size(dist.V, 1) - 1) || return false
    isRoundedPosDef(dist.V) || return false

    return true
end

function ==(x::Wishart, y::Wishart)
    is(x, y) && return true
    isApproxEqual(x.nu, y.nu) || return false
    (length(x.V)==length(y.V)) || return false
    isApproxEqual(x.V, y.V) || return false

    return true
end

# Convert Gamma -> Wishart
convert(::Type{Wishart}, d::Gamma) = Wishart(V = [1.0/(2.0*d.b)].', nu = 2.0*d.a)

# Convert Wishart -> Gamma
function convert(::Type{Gamma}, d::Wishart{1})
    Gamma(a = d.nu/2.0, b = 1.0/(2.0*d.V[1, 1]))
end

dimensions{dims}(distribution::Wishart{dims}) = (dims, dims)

dimensions{T<:Wishart}(distribution_type::Type{T}) = (distribution_type.parameters[end], distribution_type.parameters[end])

# Entropy functional
function differentialEntropy{dims}(dist::Wishart{dims})
    0.5*(dims + 1.0)*log(det(dist.V)) +
    0.5*dims*(dims + 1.0)*log(2) +
    0.25*dims*(dims - 1.0)*log(pi) +
    sum([lgamma(0.5*(dist.nu + 1.0 - i)) for i=1:dims]) -
    0.5*(dist.nu - dims - 1.0) * sum([digamma(0.5*(dist.nu + 1.0 - i)) for i=1:dims]) +
    0.5*dist.nu*dims
end
