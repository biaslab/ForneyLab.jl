export Nonlinear

"""
Description:

    Nonlinear node modeling a nonlinear relation. Updates for
    the nonlinear node are computed through the unscented transform.
    
    For more details see "On Approximate Nonlinear Gaussian Message Passing on
    Factor Graphs", Petersen et al. 2018.

    f(out, in1) = δ(out - g(in1))

Interfaces:

    1. out
    2. in1

Construction:

    Nonlinear(out, in1, g, id=:my_node)
"""
mutable struct Nonlinear <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    g_inv::Union{Function, Nothing} # Inverse of g (optional)
    dims::Tuple # Dimension of breaker message on input interface

    function Nonlinear(out, in1, g::Function; g_inv=nothing, dims=(1,), id=ForneyLab.generateId(Nonlinear))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), g, g_inv, dims)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{Nonlinear}) = "g"