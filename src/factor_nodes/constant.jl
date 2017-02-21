export Constant

"""
Description:

    A factor that clamps a variable to a constant value.

    f(x, value) = δ(x - value)

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
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        addNode!(currentGraph(), self)

        return out
    end
end

Constant(value::Any; id=generateNodeId(Constant)) = Constant(Variable(id=id), value, id=id)