############################################
# StudentsTDistribution
############################################
# Description:
#   Encodes a student's t-distribution.
#   Pamameters: m (mean), lambda (inverse scale), nu (degrees of freedom)
############################################
export StudentsTDistribution

type StudentsTDistribution <: ProbabilityDistribution
    m::Vector{Float64}      # mean
    lambda::Matrix{Float64} # inverse scale
    nu::Float64             # degrees of freedom
end

function StudentsTDistribution(; m::Union(Float64, Vector{Float64}) = [0.0],
                                 lambda::Union(Float64, Matrix{Float64}) = reshape([1.0], 1, 1),
                                 nu::Float64 = huge)
    m = (typeof(m)==Float64) ? [m] : deepcopy(m)
    lambda = (typeof(lambda)==Float64) ? fill!(Array(Float64,1,1),lambda) : deepcopy(lambda)

    return StudentsTDistribution(m, lambda, nu)
end

isProper(dist::StudentsTDistribution) = (isRoundedPosDef(dist.lambda) && (dist.nu >= tiny))

function Base.mean(dist::StudentsTDistribution)
    if isProper(dist) && dist.nu > 1
        return dist.m
    else
        return fill!(similar(dist.m), NaN)
    end
end

function Base.var(dist::StudentsTDistribution)
    if isProper(dist) && dist.nu > 2
        return dist.nu / (dist.nu - 2) * inv(dist.lambda)
    elseif isProper(dist) && dist.nu > 1
        return diagm(fill!(similar(dist.m), huge))
    else
        return fill!(similar(dist.lambda), NaN)
    end
end

format(dist::StudentsTDistribution) = "St(μ=$(format(dist.m)), λ=$(format(dist.lambda)), ν=$(format(dist.nu)))"
show(io::IO, dist::StudentsTDistribution) = println(io, format(dist))

function ==(x::StudentsTDistribution, y::StudentsTDistribution)
    return (is(x, y) || (x.m==y.m && x.lambda==y.lambda && x.nu==y.nu))
end

vague(::Type{StudentsTDistribution}) = StudentsTDistribution(m=0.0, lambda=tiny, nu=huge)
