mutable struct ScoreData{T <: String}
    id::T
    value::Real
    type::T

    ScoreData{T}(id::T, type::T, value::Real) where {T<:String} = new()
end
