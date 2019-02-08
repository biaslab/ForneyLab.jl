export GainEquality

"""
Description:

    Gain equality factor node.
    f(out,in1,in2) = δ(in2 - c'out)δ(in2 - c'in1)
    where unit vector
    c = [1, 0, ..., 0]

Interfaces:

    1. out
    2. in1
    3. in2

Construction:

    GainEquality(out, in1, in2, id=:some_id)
"""
mutable struct GainEquality <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function GainEquality(out, in1, in2, gain; id=generateId(GainEquality))
        @ensureVariables(out, in1, in2, gain)
        self = new(id, Array{Interface}(undef, 4), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)
        self.i[:gain] = self.interfaces[4] = associate!(Interface(self), gain)
        return self
    end
end

slug(::Type{GainEquality}) = "gain_eq"
