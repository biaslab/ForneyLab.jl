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
end
function StudentsTDistribution(; m::Union(Float64, Vector{Float64})=[0.0],
                                 W::Union(Float64, Matrix{Float64})=reshape([1.0], 1, 1),
                                 nu::Float64=1.0)
    m = (typeof(m)==Float64) ? [m] : deepcopy(m)
    W = (typeof(W)==Float64) ? fill!(Array(Float64,1,1),W) : deepcopy(W)

    return StudentsTDistribution(m, W, nu)
end

Base.mean(dist::StudentsTDistribution) = dist.m
Base.var(dist::StudentsTDistribution) = dist.nu / (dist.nu - 2) * inv(dist.W)

function show(io::IO, dist::StudentsTDistribution)
    println(io, typeof(dist))
    println(io, "m = $(dist.m) (location)")
    println(io, "W = $(dist.W) (precision)")
    println(io, "nu = $(dist.nu) (degrees of freedom)")
end

==(x::StudentsTDistribution, y::StudentsTDistribution) = (x.m==y.m && x.W==y.W && x.nu==y.nu)