export Sigmoid

"""
Description:
    Constrains a continuous, real-valued variable with a binary (boolean) variable.

    f(x,y) = σ(x⋅y)

Interfaces:
    1. real
    2. bin

Construction:
    Sigmoid(id=:some_id)
"""
type Sigmoid <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Sigmoid(bin::Variable, real::Variable; id=generateId(Sigmoid))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:real] = self.interfaces[1] = associate!(Interface(self), real)
        self.i[:bin] = self.interfaces[2] = associate!(Interface(self), bin)

        return self
    end
end

slug(::Type{Sigmoid}) = "σ"