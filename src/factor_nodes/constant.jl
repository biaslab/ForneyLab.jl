export Constant, constant, placeholder

"""
Description:

    A factor that clamps a variable to a constant value.

    f(x) = δ(x - value)

Interfaces:

    1. out

Construction:

    Constant(out, value, id=:some_id)
    Constant(value, id=:some_id)
"""
type Constant <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}
    value::Any

    function Constant(out::Variable, value::Any; id=generateId(Constant))
        self = new(id, Array(Interface, 1), Dict{Symbol,Interface}(), value)
        addNode!(currentGraph(), self)

        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end

"""
`constant` creates a `Variable` which is linked to a new `Constant`,
and returns this variable.

    y = constant(3.0, id=:y)
"""
function constant(value::Any; id=generateId(Constant))
    # This is basically an outer constructor for Constant,
    # but the function is not capitalized because it returns the
    # new variable linked to the Constant.
    var = Variable(id=id)
    Constant(var, value, id=id)

    return var
end

"""
`placeholder(...)` creates a `Constant` node and registers this node
as a data placeholder with the current graph.

    # Link variable y to buffer with id :y,
    # indicate that Constant will hold Float64 values.
    placeholder(y, :y, datatype=Float64)

    # Link variable y to index 3 of buffer with id :y.
    # Specify the data type by passing a default value for the Constant.
    placeholder(y, :y, index=3, default=0.0)

    # Indicate that the Constant will hold an array of size `dims`,
    # with Float64 elements.
    placeholder(X, :X, datatype=Float64, dims=(3,2))
"""
function placeholder(   var::Variable,
                        buffer_id::Symbol;
                        index::Int=0,
                        default::Any=nothing,
                        datatype::DataType=Float64,
                        dims::Tuple=())

    # Build placeholder id
    constant_id = :placeholder_*buffer_id
    if index > 0
        constant_id *= :_*index
    end

    # Define a default value
    if default != nothing
        value = default
    elseif isempty(dims)
        value = zero(datatype)
    else
        value = zeros(datatype, dims)
    end

    node = Constant(var, value, id=constant_id)
    current_graph.placeholders[node] = (buffer_id, index)

    return var
end