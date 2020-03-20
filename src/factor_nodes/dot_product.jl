export DotProduct, dot

"""
Description:

    out = in1'*in2

    in1: d-dimensional vector
    in2: d-dimensional vector
    out: scalar

           in2
           |
      in1  V   out
    ----->[⋅]----->

    f(out, in1, in2) =  δ(out - in1'*in2)

Interfaces:

    1 i[:out], 2 i[:in1], 3 i[:in2]

Construction:

    DotProduct(out, in1, in2, id=:my_node)
"""
mutable struct DotProduct <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function DotProduct(out, in1, in2; id=generateId(DotProduct))
        @ensureVariables(out, in1, in2)
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

slug(::Type{DotProduct}) = "dot"

function dot(a::Variable, in1::Variable)
    out = Variable()
    DotProduct(out, in1, a)
    return out
end

function dot(a::Variable, in1)
    @ensureVariables(in1)
    dot(a, in1)
end

function dot(a, in1::Variable)
    @ensureVariables(a)
    dot(a, in1)
end