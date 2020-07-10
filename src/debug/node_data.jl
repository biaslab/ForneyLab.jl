export NodeData

struct NodeData
    id::Union{String, Symbol}
    label::Union{String, Symbol}
    type::String
    class::String
    value::Any
end
