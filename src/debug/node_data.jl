export NodeData

mutable struct NodeData{T <: Union{String, Symbol}
    id::T
    label::T
    type::T
    class::T
    value::T

    NodeData{T}(id, label, type, class, value) where {T<:String} = new()
end
