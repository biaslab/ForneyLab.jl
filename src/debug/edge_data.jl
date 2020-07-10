export EdgeData

struct EdgeData
    id::Union{String, Symbol}
    label::String
    a::Union{String, Symbol}
    b::Union{String, Symbol}
    source::Union{String, Symbol}
    target::Union{String, Symbol}
end
