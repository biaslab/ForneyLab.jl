export EdgeData

mutable struct EdgeData{T<:Union{String, Symbol}}
    id::T
    label::T
    a::T
    b::T
    source::T
    target::T

    EdgeData{T}(id, label, a, b, source, target) where {T<:Union{String, Symbol}} = new()
end
