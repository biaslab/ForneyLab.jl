export Multiplication

import Base: *
export *

"""
Description:

    A multiplication constraint factor node

    f(x,a,z) = δ(a*x - z)

Interfaces:

    1. in
    2. gain
    3. out

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
        self.i[:in] = self.interfaces[1] = associate!(Interface(self), in1)
        self.i[:gain] = self.interfaces[2] = associate!(Interface(self), gain)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{Multiplication}) = "×"

function *(gain::Variable, in1::Variable)
    out = Variable()
    Multiplication(out, in1, gain)
    return out
end