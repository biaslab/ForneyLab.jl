export MvLogNormalDistribution

"""
Description:

    Encodes a multivariate log-normal PDF.
    Pamameters: vector m (location) and matrix S (scale).

Parameters:

    m (location vector), S (scale matrix)

Construction:

    MvLogNormalDistribution(m=zeros(3), S=eye(3))
    MvLogNormalDistribution(m=zeros(3), S=diageye(3))

Reference:

    Lognormal distributions: theory and aplications; Crow, 1988
"""
type MvLogNormalDistribution{dims} <: MultivariateProbabilityDistribution
    m::Vector{Float64} # Location
    S::AbstractMatrix{Float64} # Scale

    function MvLogNormalDistribution(m, S)
        (length(m) == size(S,1) == size(S,2)) || error("Dimensions of m and S must agree")
        return new{length(m)}(m, S)
    end
end

MvLogNormalDistribution(; m=[0.0], S=reshape([1.0],1,1)) = MvLogNormalDistribution{length(m)}(m, S)

function vague!{dims}(dist::MvLogNormalDistribution{dims})
    dist.m = zeros(dims)
    dist.S = huge*diageye(dims)
    return dist
end

vague{dims}(::Type{MvLogNormalDistribution{dims}}) = MvLogNormalDistribution(m=zeros(dims), S=huge*diageye(dims))

isProper(dist::MvLogNormalDistribution) = isRoundedPosDef(dist.S)

function Base.mean(dist::MvLogNormalDistribution)
    if isProper(dist)
        return exp(dist.m + 0.5*diag(dist.S))
    else
        return fill!(similar(dist.m), NaN)
    end
end

function Base.cov(dist::MvLogNormalDistribution)
    if isProper(dist)
        dims = size(dist, 1)
        C = zeros(dims, dims)
        for i = 1:dims
            for j = 1:dims
                C[i,j] = exp(dist.m[i] + dist.m[j] + 0.5*(S[i,i] + S[j,j]))*(exp(S[i,j]) - 1.0)
            end
        end
        return C
    else
        return fill!(similar(dist.S, NaN))
    end
end

Base.var(dist::MvLogNormalDistribution) = exp(2.0*dist.m + diag(dist.S)).*(exp(diag(dist.S)) - 1.0)

format(dist::MvLogNormalDistribution) = "logN(μ=$(format(dist.m)), Σ=$(format(dist.S)))"

show(io::IO, dist::MvLogNormalDistribution) = println(io, format(dist))

==(x::MvLogNormalDistribution, y::MvLogNormalDistribution) = (x.m==y.m && x.S==y.S)
