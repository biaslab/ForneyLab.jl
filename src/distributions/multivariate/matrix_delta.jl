export MatrixDeltaDistribution

"""
Description:

    Encodes a matrix delta distribution.
    p(X) = δ(X-M)

Parameters:

    M (location) (Number matrix)
    
Construction:
    
    MatrixDeltaDistribution(eye(3))
    MatrixDeltaDistribution(Diagonal(ones(3)))
"""
type MatrixDeltaDistribution{T, dims_n, dims_m} <: ForneyLab.MultivariateProbabilityDistribution
    M::AbstractMatrix{T}
end

MatrixDeltaDistribution{T<:Number}(M::AbstractMatrix{T}) = MatrixDeltaDistribution{T, size(M, 1), size(M, 2)}(deepcopy(M))

MatrixDeltaDistribution() = MatrixDeltaDistribution{Float64, 1, 1}(eye(1))

format(dist::MatrixDeltaDistribution) = "δ(M=$(format(dist.M)))"

show(io::IO, dist::MatrixDeltaDistribution) = println(io, format(dist))

isProper(dist::MatrixDeltaDistribution) = true

Base.mean(dist::MatrixDeltaDistribution) = dist.M

Base.mean(::Type{MatrixDeltaDistribution{Float64}}, d::MatrixDeltaDistribution) = deepcopy(d) # Definition for post-processing

sample(dist::MatrixDeltaDistribution) = dist.M

sample(::Type{MatrixDeltaDistribution{Float64}}, d::MatrixDeltaDistribution) = deepcopy(d) # Definition for post-processing

==(x::MatrixDeltaDistribution, y::MatrixDeltaDistribution) = (x.M == y.M)

convert{T<:Number}(::Type{ProbabilityDistribution}, obj::AbstractMatrix{T}) = MatrixDeltaDistribution(obj)

convert{T<:Number}(::Type{MatrixDeltaDistribution}, obj::AbstractMatrix{T}) = MatrixDeltaDistribution(obj)