export UnitEquality

"""
Description:

    Unity equality factor node.
    f(out,in1,in2) = δ(in2 - c'out)δ(in2 - c'in1)
    where unit vector
    c = [1, 0, ..., 0]

Interfaces:

    1. out
    2. in1
    3. in2

Construction:

    UnitEquality(out, in1, in2, id=:some_id)
"""
mutable struct UnitEquality <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function UnitEquality(out, in1, in2; id=generateId(UnitEquality))
        @ensureVariables(out, in1, in2)
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

slug(::Type{UnitEquality}) = "unit_eq"
