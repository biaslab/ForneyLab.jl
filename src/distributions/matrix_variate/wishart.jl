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
type Wishart{dims} <: MatrixVariate
    V::AbstractMatrix{Float64}  # Scale matrix
    nu::Float64                 # Degrees of freedom
end

Wishart(; V::AbstractMatrix{Float64} = [1.0].', nu::Float64 = 1.0) = Wishart{size(V, 1)}(deepcopy(V), nu)

Wishart() = Wishart(V = [1.0].', nu = 1.0)

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

@symmetrical function prod!{dims}(x::Wishart{dims}, y::MatrixDelta{Float64,dims,dims}, z::MatrixDelta{Float64,dims,dims}=deepcopy(y))
    # Multiplication of Wishart PDF with a delta
    (z.M == y.M) || (z.M[:] = y.M)

    return z
end

function Base.mean(dist::Wishart)
    if isProper(dist)
        return dist.nu*dist.V
    else
        return fill!(similar(dist.V), NaN)
    end
end

function Base.var(dist::Wishart)
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
