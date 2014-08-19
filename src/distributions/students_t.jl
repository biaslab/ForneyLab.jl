############################################
# StudentsTDistribution
############################################
# Description:
#   Encodes a student's t-distribution.
#   Pamameters: m (mean), W (precision), nu (degrees of freedom)
############################################
export StudentsTDistribution

type StudentsTDistribution <: ProbabilityDistribution
    m::Vector{Float64} # mean
    W::Matrix{Float64} # precision
    nu::Float64 # degrees of freedom
    StudentsTDistribution(; m=[0.0], W=reshape([1.0], 1, 1), nu=1.0) = new(m, W, nu)
end

uninformative(dist_type::Type{StudentsTDistribution}) = StudentsTDistribution(m=[0.0], W=reshape([0.001], 1, 1), nu=1000)

Base.mean(dist::StudentsTDistribution) = dist.m
Base.var(dist::StudentsTDistribution) = dist.nu / (dist.nu - 2) * inv(dist.W)

function show(io::IO, dist::StudentsTDistribution)
    println(io, typeof(dist))
    println(io, "m = $(dist.m) (location)")
    println(io, "W = $(dist.W) (precision)")
    println(io, "nu = $(dist.nu) (degrees of freedom)")
end

==(x::StudentsTDistribution, y::StudentsTDistribution) = (x.m==y.m && x.W==y.W && x.nu==y.nu)