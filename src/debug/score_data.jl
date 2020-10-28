export ScoreData

struct ScoreData
    id::Union{String, Symbol}
    value::Float64
    type::String
end
