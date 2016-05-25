export StudentsT

"""
Description:

    Encodes a student's t-distribution.

Pamameters:

    Real scalars m (mean), lambda (inverse scale), nu (degrees of freedom)

Construction:

    StudentsT(m, lambda, nu)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type StudentsT <: Univariate
    m::Float64      # mean
    lambda::Float64 # inverse scale
    nu::Float64     # degrees of freedom
end

function StudentsT(; m::Float64 = 0.0,
                                 lambda::Float64 = 1.0,
                                 nu::Float64 = huge)
    return StudentsT(m, lambda, nu)
end

function pdf(dist::StudentsT, x::Float64)
    ν = dist.nu; λ = dist.lambda; m = dist.m
    C = gamma((ν+1)/2.) / (gamma(ν/2.) * (λ/(pi*ν))^2)
    return C * (1. + λ*(x-m)^2 / ν)^(-(ν+1.)/2.)
end

isProper(dist::StudentsT) = (realmin(Float64) < abs(dist.lambda) < realmax(Float64)) && (dist.nu > realmin(Float64))

function Base.mean(dist::StudentsT)
    if isProper(dist) && dist.nu > 1
        return dist.m
    else
        return NaN
    end
end

function Base.var(dist::StudentsT)
    if isProper(dist) && dist.nu > 2
        return dist.nu / ((dist.nu - 2) * dist.lambda)
    elseif isProper(dist) && dist.nu > 1 # 1 < ν ≤ 2
        return huge
    else
        return NaN
    end
end

format(dist::StudentsT) = "St(μ=$(format(dist.m)), λ=$(format(dist.lambda)), ν=$(format(dist.nu)))"

show(io::IO, dist::StudentsT) = println(io, format(dist))

function ==(x::StudentsT, y::StudentsT)
    return (is(x, y) || (x.m==y.m && x.lambda==y.lambda && x.nu==y.nu))
end

# The choice of nu influences the robustness of the estimation.
# nu -> 0 and nu -> Inf imply higher and lower robustness respectively
# see also (Fonseca, 2008; Objective Bayesian analysis for the Student-t regression model)
# and (Lange, 1989; Robust statistical modeling using the t distribution).
# The value for nu = 4 is taken from the latter reference as a value that has worked well in practice.
function vague!(dist::StudentsT)
    dist.m = 0.0
    dist.lambda = tiny
    dist.nu = 4.0
    return dist
end

@symmetrical function prod!(x::StudentsT, y::Delta{Float64}, z::Delta{Float64}=Delta(y.m))
    # Product of log-normal PDF and Delta
    z.m = y.m

    return z
end

@symmetrical function prod!(x::StudentsT, y::Delta{Float64}, z::StudentsT)
    # Product of multivariate log-normal PDF and MvDelta, force result to be mv log-normal
    z.m = y.m
    z.lambda = huge
    z.nu = huge

    return z
end
