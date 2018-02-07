export Multiplication

import Base: *
export *

"""
Description:

    A multiplication constraint factor node

    f(out, in1, in2) = δ(out - in1*in2)

Interfaces:

    1. out
    2. in1
    3. in2

Construction:

    Multiplication(out, in1, in2, id=:some_id)
"""
mutable struct Multiplication <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Multiplication(out::Variable, in1::Variable, in2::Variable; id=generateId(Multiplication))
        self = new(id, Array{Interface}(3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

slug(::Type{Multiplication}) = "×"

function *(in1::Variable, in2::Variable)
    out = Variable()
    Multiplication(out, in1, in2)
    return out
end