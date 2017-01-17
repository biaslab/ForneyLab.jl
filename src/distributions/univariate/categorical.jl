export Categorical

"""
Description:

    Encodes a categorical (generalized Bernoulli) distribution.
    The support is {1,...,k}

Parameters:

    k > 1: number of categories
    p: probability vector of length k

Construction:

    Categorical([0.75; 0.25])

"""
type Categorical{k} <: Univariate
    p::Vector{Float64}  # probability vector
end

function Categorical(p::Vector{Float64})
    (length(p) > 1) || error("Categorical distribution should have at least two categories")
    return Categorical{length(p)}(p)
end

Categorical() = Categorical([0.5; 0.5])

function pdf(dist::Categorical, x::Int64)
    return (1 <= x <= length(dist.p)) ? dist.p[x] : 0.0
end

function vague!{k}(dist::Categorical{k})
    fill!(dist.p, 1./k)
    return dist
end

vague{k}(::Type{Categorical{k}}) = Categorical((1./k)*ones(k))

format{k}(dist::Categorical{k}) = "Cat{$(k)}($(format(dist.p)))"

show(io::IO, dist::Categorical) = println(io, format(dist))

isProper(dist::Categorical) = (abs(sum(dist.p)-1.) < 1e-6)

unsafeMean(dist::Categorical) = deepcopy(dist.p)

function prod!{k}(x::Categorical{k}, y::Categorical{k}, z::Categorical{k}=vague(Categorical{k}))
    # Multiplication of 2 categorical PMFs: p(z) = p(x) * p(y)
    z.p[:] = x.p .* y.p
    norm = sum(z.p)
    (norm > 0.) || error("Product of categoricals cannot be normalized")
    z.p *= 1./norm

    return z
end

@symmetrical function prod!{k}(x::Categorical{k}, y::Delta{Int64}, z::Delta{Int64}=Delta(y.m))
    # Product of categorical PMF and Delta
    (1 <= y.m <= k) || throw(DomainError())
    (x.p[y.m] > 0.) || error("Invalid product of Categorical and Delta: result cannot be normalized")
    z.m = y.m

    return z
end

@symmetrical function prod!{k}(x::Categorical{k}, y::Delta{Int64}, z::Categorical{k})
    # Product of categorical PMF and Delta, force result to be Categorical
    (1 <= y.m <= k) || throw(DomainError())
    (x.p[y.m] > 0.) || error("Invalid product of Categorical and Delta: result cannot be normalized")
    z.p[:] = 0.
    z.p[y.m] = 1.

    return z
end

function ==(x::Categorical, y::Categorical)
    is(x, y) && return true
    if (typeof(x)==typeof(y)) && isApproxEqual(x.p, y.p)
        return true
    else
        return false
    end
end