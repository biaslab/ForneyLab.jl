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

function constant(value::Any; id=generateId(Constant))
    # This is basically an outer constructor for Constant,
    # but it returns the corresponding Variable object rather
    # than the Constant object.
    # Therefore, the function name is not capitalized.
    var = Variable(id=id)
    Constant(var, value, id=id)

    return var
end

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