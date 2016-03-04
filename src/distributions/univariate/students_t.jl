export StudentsTDistribution

"""
Description:

    Encodes a student's t-distribution.

Pamameters:

    Real scalars m (mean), lambda (inverse scale), nu (degrees of freedom)

Construction:

    StudentsTDistribution(m, lambda, nu)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
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

# The choice of nu influences the robustness of the estimation.
# nu -> 0 and nu -> Inf imply higher and lower robustness respectively
# see also (Fonseca, 2008; Objective Bayesian analysis for the Student-t regression model)
# and (Lange, 1989; Robust statistical modeling using the t distribution).
# The value for nu = 4 is taken from the latter reference as a value that has worked well in practice.
function vague!(dist::StudentsTDistribution)
    dist.m = 0.0
    dist.lambda = tiny
    dist.nu = 4.0
    return dist
end
