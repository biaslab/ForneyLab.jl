export Multiplication

import Base: *
export *

"""
Description:

    A multiplication constraint factor node

    f(out, in, gain) = δ(gain*in - out)

Interfaces:

    1. out
    2. in
    3. gain

Construction:

    Multiplication(out, in, gain, id=:some_id)
"""
type Multiplication <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Multiplication(out::Variable, in1::Variable, gain::Variable; id=generateId(Multiplication))
        self = new(id, Array(Interface, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:gain] = self.interfaces[3] = associate!(Interface(self), gain)

        return self
    end
end

slug(::Type{Multiplication}) = "×"

function *(gain::Variable, in1::Variable)
    out = Variable()
    Multiplication(out, in1, gain)
    return out
end