export MatrixDelta

"""
Description:

    Encodes a matrix delta distribution.
    p(X) = δ(X-M)

Parameters:

    M (location) (Number matrix)

Construction:

    MatrixDelta(eye(3))
    MatrixDelta(diageye(3))
"""
type MatrixDelta{T, dims_n, dims_m} <: MatrixVariate
    M::AbstractMatrix{T}
end

MatrixDelta{T<:Number}(M::AbstractMatrix{T}) = MatrixDelta{T, size(M, 1), size(M, 2)}(deepcopy(M))

MatrixDelta() = MatrixDelta{Float64, 1, 1}(eye(1))

pdf(dist::MatrixDelta, x::AbstractMatrix) = ((dist.M==x) ? 1.0 : 0.0)

format(dist::MatrixDelta) = "δ(M=$(format(dist.M)))"

show(io::IO, dist::MatrixDelta) = println(io, format(dist))

function prod!{T,dims_n,dims_m}(x::MatrixDelta{T,dims_n,dims_m}, y::MatrixDelta{T,dims_n,dims_m}, z::MatrixDelta{T,dims_n,dims_m}=deepcopy(y))
    (x.M == y.M) || throw(DomainError())
    (z.M == y.M) || (z.M = deepcopy(y.M))

    return z
end

isProper(dist::MatrixDelta) = true

Base.mean(dist::MatrixDelta) = dist.M

Base.mean(::Type{MatrixDelta{Float64}}, d::MatrixDelta) = deepcopy(d) # Definition for post-processing

sample(dist::MatrixDelta) = dist.M

sample(::Type{MatrixDelta{Float64}}, d::MatrixDelta) = deepcopy(d) # Definition for post-processing

==(x::MatrixDelta, y::MatrixDelta) = (x.M == y.M)

convert{T<:Number}(::Type{ProbabilityDistribution}, obj::AbstractMatrix{T}) = MatrixDelta(obj)

convert{T<:Number}(::Type{MatrixDelta}, obj::AbstractMatrix{T}) = MatrixDelta(obj)

dimensions{T, dims_n, dims_m}(distribution::MatrixDelta{T,dims_n,dims_m}) = (dims_n, dims_m)

dimensions{T<:MatrixDelta}(distribution_type::Type{T}) = (distribution_type.parameters[end-1], distribution_type.parameters[end])
