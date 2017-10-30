export Equality

"""
Description:

    An equality constraint factor node

    f([1],[2],[3]) = δ([1] - [2]) δ([1] - [3])

Interfaces:

    1, 2, 3

Construction:

    Equality(id=:some_id)

    The interfaces of an Equality node have to be connected manually.
"""
type Equality <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function Equality(; id=generateId(Equality))
        self = new(id, Array(Interface, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:3
            self.i[idx] = self.interfaces[idx] = Interface(self)
        end

        return self
    end
end

slug(::Type{Equality}) = "="