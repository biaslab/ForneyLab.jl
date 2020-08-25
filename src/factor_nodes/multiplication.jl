export Multiplication, *

"""
Description:

    For continuous random variables, the multiplication node acts
    as a (matrix) multiplication constraint, with node function

    f(out, in1, a) = δ(out - a*in1)

Interfaces:

    1. out
    2. in1
    3. a

Construction:

    Multiplication(out, in1, a, id=:some_id)
"""
mutable struct Multiplication <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Multiplication(out, in1, a; id=generateId(Multiplication))
        @ensureVariables(out, in1, a)
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:a] = self.interfaces[3] = associate!(Interface(self), a)

        return self
    end
end

slug(::Type{Multiplication}) = "×"

function *(a::Variable, in1::Variable)
    out = Variable()
    Multiplication(out, in1, a)
    return out
end

function *(in1::Variable, a)
    @ensureVariables(a)
    *(a, in1)
end

function *(a, in1::Variable)
    @ensureVariables(a)
    *(a, in1)
end