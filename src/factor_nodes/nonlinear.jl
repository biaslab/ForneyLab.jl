export Nonlinear

"""
Description:

    Nonlinear node modeling a nonlinear relation. Updates for 
    the nonlinear node are computed through local linearization.

    f(out, in1) = δ(out - g(in1))

Interfaces:

    1. out
    2. in1

Construction:

    Nonlinear(out, in1, g::Function, J_g::Function, id=:my_node)
"""
mutable struct Nonlinear <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    J_g::Function # Jacobi matrix of g, as a function of the input vector; in the 1-d case this reduces to the first derivative of g

    function Nonlinear(out::Variable, in1::Variable, g::Function, J_g::Function; id=ForneyLab.generateId(Nonlinear))
        self = new(id, Vector{Interface}(2), Dict{Symbol,Interface}(), g, J_g)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{Nonlinear}) = "g"