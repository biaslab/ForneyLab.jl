############################################
# StudentsTDistribution
############################################
# Description:
#   Encodes a student's t-distribution.
#   Pamameters: m (mean), lambda (inverse scale), nu (degrees of freedom)
############################################
export StudentsTDistribution

type StudentsTDistribution <: UnivariateProbabilityDistribution
    m::Float64      # mean
    lambda::Float64 # inverse scale
    nu::Float64     # degrees of freedom
end

function StudentsTDistribution(; m::Float64 = 0.0,
                                 lambda::Float64 = 1.0,
                                 nu::Float64 = huge)
    return StudentsTDistribution(m, lambda, nu)
end

isProper(dist::StudentsTDistribution) = (realmin(Float64) < abs(dist.lambda) < realmax(Float64)) && (dist.nu > realmin(Float64))

function Base.mean(dist::StudentsTDistribution)
    if isProper(dist) && dist.nu > 1
        return dist.m
    else
        return NaN
    end
end

function Base.var(dist::StudentsTDistribution)
    if isProper(dist) && dist.nu > 2
        return dist.nu / ((dist.nu - 2) * dist.lambda)
    elseif isProper(dist) && dist.nu > 1 # 1 < ν ≤ 2
        return huge
    else
        return NaN
    end
end

format(dist::StudentsTDistribution) = "St(μ=$(format(dist.m)), λ=$(format(dist.lambda)), ν=$(format(dist.nu)))"

show(io::IO, dist::StudentsTDistribution) = println(io, format(dist))

function ==(x::StudentsTDistribution, y::StudentsTDistribution)
    return (is(x, y) || (x.m==y.m && x.lambda==y.lambda && x.nu==y.nu))
end

vague(::Type{StudentsTDistribution}) = StudentsTDistribution(m=0.0, lambda=tiny, nu=huge)