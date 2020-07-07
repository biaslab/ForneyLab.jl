mutable struct NodeData{T <: String}
    id::T
    label::T
    type::T
    class::T
    value::T

    NodeData{T}(id, label, type, class, value) where {T<:String} = new()
end
