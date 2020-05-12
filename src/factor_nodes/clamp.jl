import Base: copy
export Clamp, constant, placeholder, @ensureVariables

"""
Description:

    A factor that clamps a variable to a constant value.

    f(out) = δ(out - value)

Interfaces:

    1. out

Construction:

    Clamp(out, value, id=:some_id)
    Clamp(value, id=:some_id)
"""
mutable struct Clamp{T<:VariateType} <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}
    value::Any

    # Fields for algorithm assembly
    dist_or_msg::Type # Define whether the clamp outputs a Message or ProbabilityDistribution
    buffer_id::Symbol # Specify the buffer id for an attached placeholder
    buffer_index::Int64 # Specify the buffer index when the placeholder attaches a vector

    Clamp(V::Type) = new{V}() # Instantiate an empty clamp without adding it to the graph

    function Clamp(out::Variable, value::Any; id=generateId(Clamp{variateType(value)}))
        self = new{variateType(value)}(id, Array{Interface}(undef, 1), Dict{Symbol,Interface}(), value)
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end

"""
Make a shallow copy of a Clamp without adding it to the graph.
The copied Clamp can be used as altered inboud during algorithm
assembly, without altering the original source Clamp.
"""
function copy(src::Clamp{V}) where V<:VariateType
    nd = Clamp(V)
    for field in fieldnames(Clamp)
        isdefined(src, field) && setfield!(nd, field, getfield(src, field))
    end

    return nd
end

variateType(value::Number) = Univariate
variateType(value::Vector) = Multivariate
variateType(value::AbstractMatrix) = MatrixVariate

"""
`constant` creates a `Variable` which is linked to a new `Clamp`,
and returns this variable.

    y = constant(3.0, id=:y)
"""
function constant(value::Any; id=generateId(Clamp{variateType(value)}))
    # This is basically an outer constructor for Clamp,
    # but the function returns the new variable linked to the Clamp.
    var = Variable(id=id)
    Clamp(var, value, id=id)

    return var
end

"""
`@ensureVariables(...)` casts all non-`Variable` arguments to `Variable` through constant(arg).
"""
macro ensureVariables(args...)
    checks = map(arg -> :(isa($arg, Variable) || ($arg = constant($arg))), args)
    expr   = reduce((current, item) -> :($current; $item), checks)
    return esc(expr)
end

"""
`placeholder(...)` creates a `Clamp` node and registers this node
as a data placeholder with the current graph.

    # Link variable y to buffer with id :y,
    # indicate that Clamp will hold Float64 values.
    placeholder(y, :y, datatype=Float64)

    # Link variable y to index 3 of buffer with id :y.
    # Specify the data type by passing a default value for the Clamp.
    placeholder(y, :y, index=3, default=0.0)

    # Indicate that the Clamp will hold an array of size `dims`,
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
    constant_id = :placeholder_*var.id

    # Define a default value
    if default != nothing
        value = default
    elseif isempty(dims)
        value = zero(datatype)
    else
        value = zeros(datatype, dims)
    end

    node = Clamp(var, value, id=constant_id)
    current_graph.placeholders[node] = (buffer_id, index)

    return var
end

function placeholder(   buffer_id::Symbol;
                        index::Int=0,
                        default::Any=nothing,
                        datatype::DataType=Float64,
                        dims::Tuple=(),
                        var_id=buffer_id)

    var = Variable(id=var_id)
    return placeholder(var, buffer_id, index=index, default=default, datatype=datatype, dims=dims)
end
