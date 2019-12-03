export Addition, +, -

"""
Description:

    An addition constraint factor node.

    f(out,in1,in2) = δ(in1 + in2 - out)

Interfaces:

    1. out
    2. in1
    3. in2

Construction:

    Addition(out, in1, in2, id=:some_id)
"""
mutable struct Addition <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Addition(out, in1, in2; id=generateId(Addition))
        @ensureVariables(out, in1, in2)
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

slug(::Type{Addition}) = "+"

function +(in1::Variable, in2::Variable)
    out = Variable()
    Addition(out, in1, in2)
    return out
end

function +(in1::Variable, in2)
    @ensureVariables(in2)
    in1 + in2
end

+(in1, in2::Variable) = +(in2, in1)


"""
    -(in1::Variable, in2::Variable)

A subtraction constraint based on the addition factor node.

"""
function -(in1::Variable, in2::Variable)
    out = Variable()
    Addition(in1, out, in2)
    return out
end

function -(in1::Variable, in2)
    @ensureVariables(in2)
    in1 - in2
end

function -(in1, in2::Variable)
    @ensureVariables(in1)
    in1 - in2
end
