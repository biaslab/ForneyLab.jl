export Gamma

"""
Description:

    A gamma node with shape-rate parameterization:

    f(x,a,b) = Gam(x|a,b)

Interfaces:

    1. shape
    2. rate
    3. out

Construction:

    Gamma(out, shape, rate, id=:some_id)
"""
type Gamma <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Gamma(out::Variable, shape::Variable, rate::Variable; id=generateId(Gamma))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:shape] = self.interfaces[1] = associate!(Interface(self), shape)
        self.i[:rate] = self.interfaces[2] = associate!(Interface(self), rate)
        self.i[:out] = self.interfaces[3] = associate!(Interface(self), out)

        return self
    end
end

slug(::Type{Gamma}) = "Gam"