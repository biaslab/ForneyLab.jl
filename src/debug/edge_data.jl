mutable struct EdgeData{T<:String}
    id::T
    label::T
    a::T
    b::T
    source::T
    target::T

    EdgeData{T}(id, label, a, b, source, target) where {T<:String} = new()
end
