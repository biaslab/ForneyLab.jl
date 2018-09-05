export Equality, equal

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
mutable struct Equality <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function Equality(; id=generateId(Equality))
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:3
            self.i[idx] = self.interfaces[idx] = Interface(self)
        end

        return self
    end

    function Equality(out, in1, in2; id=generateId(Equality))
        # Explicit node constructor for creating an equality node where all connected
        # variables are considered separately
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[1] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[2] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[3] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

slug(::Type{Equality}) = "="

function equal(in1::Variable, in2::Variable)
    out = Variable()
    Equality(out, in1, in2)
    return out
end