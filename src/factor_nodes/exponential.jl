export Exponential, exp

"""
Description:

    Maps a location to a scale parameter by exponentiation

    f(out,in1) = δ(out - exp(in1))

Interfaces:

    1. out
    2. in1

Construction:

    Exponential(out, in1, id=:some_id)
"""
mutable struct Exponential <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Exponential(out, in1; id=generateId(Exponential))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{Exponential}) = "exp"

function exp(in1::Variable)
    out = Variable()
    Exponential(out, in1)
    return out
end