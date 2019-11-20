export Nonlinear, NonlinearUT, NonlinearPT

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
    out::Variable
    in1::Variable
    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d

    function Nonlinear(out, in1, g::Function; id=ForneyLab.generateId(Nonlinear))
        @ensureVariables(out, in1)
        self = new(id, out, in1, g)
        return self
    end
end

slug(::Type{Nonlinear}) = "g"

mutable struct NonlinearUT <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    g_inv::Union{Function, Nothing} # Inverse of g (optional)
    alpha::Union{Float64, Nothing} # Spread parameter for unscented transform
    dims::Tuple # Dimension of breaker message on input interface

    function NonlinearUT(out, in1, g::Function; g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(Nonlinear))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), g, g_inv, alpha, dims)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end

    function NonlinearUT(node::Nonlinear; g_inv=nothing, alpha=nothing, dims=())
        self = new(node.id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), node.g, g_inv, alpha, dims)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), node.out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), node.in1)

        return self
    end
end

slug(::Type{NonlinearUT}) = "g"

"""
Description:

    Deterministic functions that cause non-conjugacy can be incorporated into
    the existing factor graph with a Nonconjugate node.

Interfaces:

    1. out
    2. in1

Construction:

    NonlinearPT(out, in1, g::Function)
"""
mutable struct NonlinearPT <: DeltaFactor

    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol, Interface}

    g::Function # Deterministic function of output variable

    function NonlinearPT(out, in1, g::Function, id=ForneyLab.generateId(NonlinearPT))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}(), g)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end

    function NonlinearPT(node::Nonlinear)
        self = new(node.id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), node.g)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), node.out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), node.in1)

        return self
    end
end

slug(::Type{NonlinearPT}) = "g"
