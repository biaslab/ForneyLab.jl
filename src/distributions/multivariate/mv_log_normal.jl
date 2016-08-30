export MvLogNormal

"""
Description:

    Encodes a multivariate log-normal PDF.
    Pamameters: vector m (location) and matrix S (scale).

Parameters:

    m (location vector), S (scale matrix)

Construction:

    MvLogNormal(m=zeros(3), S=eye(3))
    MvLogNormal(m=zeros(3), S=diageye(3))

Reference:

    Lognormal distributions: theory and aplications; Crow, 1988
"""
type MvLogNormal{dims} <: Multivariate
    m::Vector{Float64} # Location
    S::AbstractMatrix{Float64} # Scale

    function MvLogNormal(m, S)
        (length(m) == size(S,1) == size(S,2)) || error("Dimensions of m and S must agree")
        return new{length(m)}(m, S)
    end
end

MvLogNormal(; m=[0.0], S=reshape([1.0],1,1)) = MvLogNormal{length(m)}(m, S)

function pdf{dims}(dist::MvLogNormal{dims}, x::Vector{Float64})
    (length(x) == dims) || return 0.0
    all(x .>= 0.0) || return 0.0
    error("TODO: pdf(::MvLogNormal, x) is not yet implemented")
end

function vague!{dims}(dist::MvLogNormal{dims})
    dist.m = zeros(dims)
    dist.S = huge*diageye(dims)
    return dist
end

vague{dims}(::Type{MvLogNormal{dims}}) = MvLogNormal(m=zeros(dims), S=huge*diageye(dims))

isProper(dist::MvLogNormal) = isRoundedPosDef(dist.S)

unsafeMean(dist::MvLogNormal) = exp(dist.m + 0.5*diag(dist.S))

function unsafeCov(dist::MvLogNormal)
    dims = dimensions(dist)
    C = zeros(dims, dims)
    for i = 1:dims
        for j = 1:dims
            C[i,j] = exp(dist.m[i] + dist.m[j] + 0.5*(S[i,i] + S[j,j]))*(exp(S[i,j]) - 1.0)
        end
    end
    return C
end

unsafeVar(dist::MvLogNormal) = exp(2.0*dist.m + diag(dist.S)).*(exp(diag(dist.S)) - 1.0)

format(dist::MvLogNormal) = "logN(μ=$(format(dist.m)), Σ=$(format(dist.S)))"

show(io::IO, dist::MvLogNormal) = println(io, format(dist))

==(x::MvLogNormal, y::MvLogNormal) = (x.m==y.m && x.S==y.S)

@symmetrical function prod!{dims}(x::MvLogNormal{dims}, y::MvDelta{Float64,dims}, z::MvDelta{Float64,dims}=MvDelta(y.m))
    # Product of multivariate log-normal PDF and MvDelta
    all(y.m .>= 0.) || throw(DomainError())
    (z.m == y.m) || (z.m[:] = y.m)

    return z
end

@symmetrical function prod!{dims}(x::MvLogNormal{dims}, y::MvDelta{Float64,dims}, z::MvLogNormal{dims})
    # Product of multivariate log-normal PDF and MvDelta, force result to be mv log-normal
    all(y.m .>= 0.) || throw(DomainError())
    z.m[:] = log(y.m)
    z.S = tiny.*eye(dims)

    return z
end
