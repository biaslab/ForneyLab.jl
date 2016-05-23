export Dirichlet

"""
Description:

    Encodes a Dirichlet distribution.

Parameters:

    alpha: vector of positive reals of length k

Construction:

    Dirichlet([1.0; 3.0; 2.0])
"""
type Dirichlet{dims} <: Multivariate
    alpha::Vector{Float64}
end

Dirichlet(alpha::Vector{Float64}) = Dirichlet{length(alpha)}(alpha)

Dirichlet() = Dirichlet{2}(ones(2))

function pdf(dist::Dirichlet, x::Vector{Float64})
    (length(x) == length(dist.alpha)) || return 0.0
    all(0.0 .<= x .<= 1.0) || return 0.0
    isApproxEqual(sum(x), 1.0) || return 0.0
    α = dist.alpha
    B_inv = gamma(sum(α)) / prod(gamma(α))

    return B_inv * prod(x.^(α-1.0))
end

format(dist::Dirichlet) = "Dir(α=$(format(dist.alpha)))"

show(io::IO, dist::Dirichlet) = println(io, format(dist))

isProper(dist::Dirichlet) = all(dist.alpha .> 0.)

function vague!{dims}(dist::Dirichlet{dims})
    fill!(dist.alpha, tiny)

    return dist
end

vague{dims}(::Type{Dirichlet{dims}}) = Dirichlet{dims}(zeros(dims)+tiny)

function ==(x::Dirichlet, y::Dirichlet)
    is(x,y) && return true
    if (typeof(x)==typeof(y)) && isApproxEqual(x.alpha, y.alpha)
        return true
    else
        return false
    end
end

function prod!{dims}(x::Dirichlet{dims}, y::Dirichlet{dims}, z::Dirichlet{dims}=deepcopy(y))
    # Multiplication of 2 Dirichlet PDFs: p(z) = p(x) * p(y)
    z.alpha = x.alpha + y.alpha - 1.0

    return z
end

@symmetrical function prod!{dims}(x::Dirichlet{dims}, y::MvDelta{Float64,dims}, z::MvDelta{Float64,dims}=deepcopy(y))
    # Multiplication of a Dirichlet PDF with a MvDelta
    all(0.0 .<= y.m .<= 1.0) || throw(DomainError())
    isApproxEqual(sum(y.m), 1.0) || throw(DomainError())
    z.m[:] = y.m

    return z
end

@symmetrical function prod!{dims}(x::Dirichlet{dims}, y::MvDelta{Float64,dims}, z::Dirichlet{dims})
    # Multiplication of a Dirichlet PDF with a MvDelta, force result to be Dirichlet
    error("TODO")
end