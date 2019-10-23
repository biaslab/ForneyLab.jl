export Nonconjugate

"""
Description:

    Deterministic functions that cause non-conjugacy can be incorporated into
    the existing factor graph with Nonconjugate node. Primitive implementation requires input
    to be Gaussian distributed random variable and approximates the posterior marginal by
    Gauss-Hermite quadrature/cubature technique. The node computes a pdf function by change of random variables
    and autodiff, then transmits the pdf message through in1 interface. When it is needed out interface carries
    an abstract distribution message which is not a parametric distribution but supplies the
    required sufficient statistics to carry out VMP and free energy computation.

Interfaces:

    1. out
    2. in1

Construction:

    Nonconjugate(out, in1, g::Function)
"""
mutable struct Nonconjugate <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol, Interface}

    g::Function # Deterministic function of output variable

    function Nonconjugate(out, in1, g::Function, id=ForneyLab.generateId(Nonconjugate))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}(), g)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{Nonconjugate}) = "Nonconjugate"
