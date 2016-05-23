export PartitionedDistribution

"""
Description:

    Encodes the product of multiple independent distributions of the same type.
    `p(x1,x2,x3) = p(x1)p(x2)p(x3)`

Parameters:

    dtype<:ProbabilityDistribution: type of the factors
    n_factors: number of factors

Fields:

    factors::Vector{ProbabilityDistribution}

Construction:

    PartitionedDistribution([Gaussian(), Gaussian()])
"""
type PartitionedDistribution{dtype<:ProbabilityDistribution,n_factors} <: Multivariate
    factors::Vector{dtype}

    function PartitionedDistribution(factors::Vector{dtype})
        (length(factors) > 1) || error("PartitionedDistribution should contain at least 2 factors")
        for i=2:length(factors)
            (typeof(factors[i]) == dtype) || error("All factors in PartitionedDistribution should have the same type")
        end

        return new{dtype,length(factors)}(factors)
    end
end

PartitionedDistribution{T<:ProbabilityDistribution}(factors::Vector{T}) = PartitionedDistribution{typeof(factors[1]), length(factors)}(factors)

PartitionedDistribution() = PartitionedDistribution([Gaussian(), Gaussian()])

function pdf(dist::PartitionedDistribution, x::Vector)
    (length(x) == dimensions(dist)) || throw(DimensionMismatch())
    d = dimensions(dist.factors[1])
    p = 1.0
    factor_idx = 1
    factor_dim = 1
    for i=1:length(x)
        p *= pdf(dist.factors[factor_idx], x[i])
        if factor_dim == d
            factor_dim = 1
            factor_idx += 1
        else
            factor_dim += 1
        end
    end

    return p
end

function show(io::IO, dist::PartitionedDistribution)
    println(io, typeof(dist))
    println(io, dist.factors)
end

function vague!{dtype,n_factors}(dist::PartitionedDistribution{dtype,n_factors})
    map(vague!, dist.factors)

    return dist
end

vague{dtype,n_factors}(::Type{PartitionedDistribution{dtype,n_factors}}) = PartitionedDistribution(dtype[vague(dtype) for i=1:n_factors])

Base.mean(dist::PartitionedDistribution) = vcat(map(mean, dist.factors)...)

sample(dist::PartitionedDistribution) = vcat(map(sample, dist.factors)...)

isProper(dist::PartitionedDistribution) = all(map(isProper, dist.factors))

dimensions{dtype,n_factors}(distribution::PartitionedDistribution{dtype,n_factors}) = dimensions(dtype)*n_factors

dimensions{dtype,n_factors}(distribution_type::Type{PartitionedDistribution{dtype,n_factors}}) = dimensions(dtype)*n_factors

function prod!{Tx,Ty,Tz,n_factors}(x::PartitionedDistribution{Tx,n_factors}, y::PartitionedDistribution{Ty,n_factors}, z::PartitionedDistribution{Tz,n_factors}=vague(PartitionedDistribution{Tx,n_factors}))
    # Multiplication of 2 partitioned PDFs: p(z1)...p(zn) = p(x1)*p(y1)...p(xn)*p(yn)
    for f=1:n_factors
        prod!(x.factors[f], y.factors[f], z.factors[f])
    end

    return z
end

@symmetrical function prod!{dtype,n_factors}(x::PartitionedDistribution{dtype,n_factors}, y::MvDelta, z::MvDelta=MvDelta(zeros(2)))
    # Multiplication of a partitioned PDF with a delta
    dims = dimensions(x)
    factor_dims = dimensions(x.factors[1])
    (dims == dimensions(y)) || throw(DimensionMismatch())
    (length(z.m) == dims) || (z.m = zeros(dims))
    for i=1:n_factors
        if dtype <: Multivariate
            idx = 1 + (i-1)*factor_dims
            factor_range = idx:idx+factor_dims-1
            y_factor = MvDelta(y.m[factor_range])
            z.m[factor_range] = (x.factors[i] * y_factor).m
        else
            y_factor = Delta(y.m[i])
            z.m[i] = (x.factors[i] * y_factor).m
        end
    end

    return z
end

function ==(x::PartitionedDistribution, y::PartitionedDistribution)
    (length(x.factors)==length(x.factors)) || return false
    return all(i->(x.factors[i]==y.factors[i]), collect(1:length(x.factors)))
end
