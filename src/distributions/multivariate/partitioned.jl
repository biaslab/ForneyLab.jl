export Partitioned

"""
Description:

    Encodes the product of multiple independent distributions of the same type.
    `p(x1,x2,x3) = p(x1)p(x2)p(x3)`

Parameters:

    dtype<:ProbabilityDistribution: uniform type of the factors
    dtype<:Union{ProbabilityDistribution}: type for factors of distinct distributions
    n_factors: number of factors

Fields:

    factors::Vector{ProbabilityDistribution}

Construction:

    Partitioned([Gaussian(), Gaussian()])
"""
type Partitioned{dtype, n_factors} <: Multivariate
    factors::Vector{ProbabilityDistribution}
end

function Partitioned{T<:ProbabilityDistribution}(factors::Vector{T})
    # Construct distribution type parametrization
    unique_partition_types = unique(map(typeof, factors))
    if length(unique_partition_types) == 1 # factor distribution types are uniform
        dtype = unique_partition_types[1]
    else # factors consist of varying distribution types
        dtype = Union{unique_partition_types...}
    end

    return Partitioned{dtype, length(factors)}(factors)
end

Partitioned() = Partitioned([Gaussian()])

function pdf{dtype,n_factors}(dist::Partitioned{dtype,n_factors}, x::Vector)
    (length(x) == dimensions(dist)) || throw(DimensionMismatch())

    p = 1.0
    x_start_idx = 1 # starting position for x
    for factor_idx = 1 : n_factors # iterate over factors
        factor_dist = dist.factors[factor_idx] # get the factor distribution
        factor_dims = dimensions(factor_dist) # get the factor distribution dimensionality

         # get the corresponding x indices/index
        if factor_dims > 1
            x_idc = x_start_idx : x_start_idx+factor_dims-1
        else
            x_idc = x_start_idx
        end

        p *= pdf(factor_dist, x[x_idc]) # evaluate the factor distribution at the appropriate values obtained from x
        x_start_idx += factor_dims
    end

    return p
end

function show(io::IO, dist::Partitioned)
    println(io, typeof(dist))
    println(io, dist.factors)
end

function vague!{dtype,n_factors}(dist::Partitioned{dtype,n_factors})
    map(vague!, dist.factors)

    return dist
end

vague{dtype<:ProbabilityDistribution,n_factors}(::Type{Partitioned{dtype,n_factors}}) = Partitioned(dtype[vague(dtype) for i=1:n_factors]) # Note: this function is only implemented for a uniform partition

m(dist::Partitioned) = vcat(map(m, dist.factors)...)

sample(dist::Partitioned) = vcat(map(sample, dist.factors)...)

isProper(dist::Partitioned) = all(map(isProper, dist.factors))

dimensions{dtype,n_factors}(distribution::Partitioned{dtype,n_factors}) = sum(map(dimensions, distribution.factors))

dimensions{dtype<:ProbabilityDistribution, n_factors}(distribution_type::Type{Partitioned{dtype,n_factors}}) = dimensions(dtype)*n_factors # Note: this function is only implemented for a uniform partition

function prod!{Tx,Ty,Tz,n_factors}(x::Partitioned{Tx,n_factors}, y::Partitioned{Ty,n_factors}, z::Partitioned{Tz,n_factors}=vague(Partitioned{Tx,n_factors}))
    # Multiplication of 2 partitioned PDFs: p(z1)...p(zn) = p(x1)*p(y1)...p(xn)*p(yn)
    for f=1:n_factors
        prod!(x.factors[f], y.factors[f], z.factors[f])
    end

    return z
end

@symmetrical function prod!{dtype,n_factors}(x::Partitioned{dtype,n_factors}, y::MvDelta, z::MvDelta=MvDelta(zeros(dimensions(y))))
    # Multiplication of a partitioned PDF with a delta
    (dimensions(x) == dimensions(y) == dimensions(z)) || throw(DimensionMismatch())

    y_start_idx = 1
    for factor_idx = 1 : n_factors
        factor_dist = x.factors[factor_idx]
        factor_dims = dimensions(factor_dist)

        # Select whether to multiply with a multi- or univariate delta
        if factor_dims > 1
            y_idc = y_start_idx : y_start_idx+factor_dims-1
            y_factor_dist = MvDelta(y.m[y_idc])
        else
            y_idc = y_start_idx
            y_factor_dist = Delta(y.m[y_idc])
        end

        z.m[y_idc] = (factor_dist * y_factor_dist).m
        y_start_idx += factor_dims
    end

    return z
end

function ==(x::Partitioned, y::Partitioned)
    (length(x.factors)==length(x.factors)) || return false
    return all(i->(x.factors[i]==y.factors[i]), collect(1:length(x.factors)))
end
