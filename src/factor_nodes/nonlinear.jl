export Nonlinear, Unscented, Sampling

abstract type ApproximationMethod end
abstract type Unscented <: ApproximationMethod end
abstract type Sampling <: ApproximationMethod end

"""
Description:

    Nonlinear node modeling a nonlinear relation. Updates for
    the nonlinear node are computed through the unscented transform (by default) or using importance sampling.

    For more details see "On Approximate Nonlinear Gaussian Message Passing on
    Factor Graphs", Petersen et al. 2018.

    f(out, in1) = δ(out - g(in1))

Interfaces:

    1. out
    2. in1

Construction:

    Nonlinear(out, in1; g=g, id=:my_node)
    Nonlinear(out, in1; g=g, g_inv=g_inv, id=:my_node)
    Nonlinear(out, in1, in2, ...; g=g, id=:my_node)
    Nonlinear(out, in1, in2, ...; g=g, g_inv=(g_inv_in1, g_inv_in2, ...), id=:my_node)
    Nonlinear(out, in1, in2, ...; g=g, g_inv=(g_inv_in1, nothing, ...), id=:my_node)
"""
mutable struct Nonlinear{T<:ApproximationMethod} <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output as a function of the inputs
    g_inv::Union{Function, Nothing, Vector} # Inverse of g with respect to individual inbounds (optional)
    alpha::Union{Float64, Nothing} # Spread parameter for unscented transform
    dims::Union{Tuple, Vector} # Dimension of breaker message(s) on input interface(s)
    status::Dict{Symbol, Union{Bool, Message}} #Keeps the status of node to ensure that input variables are updated simultaneously
    n_samples::Int # Number of samples for sampling

    function Nonlinear{Unscented}(out, args::Vararg; g::Function, g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(Nonlinear{Unscented}))
        @ensureVariables(out)
        n_args = length(args)
        for i=1:n_args
            @ensureVariables(args[i])
        end
        self = new(id, Vector{Interface}(undef, n_args+1), Dict{Symbol,Interface}(), g, g_inv, alpha, dims, Dict{Symbol,Bool}(), 0)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        for k = 1:n_args
            self.i[:in*k] = self.interfaces[k+1] = associate!(Interface(self), args[k])
        end

        return self
    end

    function Nonlinear{Sampling}(out, in1, g::Function; n_samples=1000, id=ForneyLab.generateId(Nonlinear{Sampling}))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), g, nothing, nothing, (), Dict{Symbol,Bool}(), n_samples)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end

    function Nonlinear{Sampling}(out, in1, in2, g::Function; n_samples=1000, id=ForneyLab.generateId(Nonlinear{Sampling}))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 3), Dict{Symbol,Interface}(), g, nothing, nothing, (), Dict(:updated=>false), n_samples)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

function Nonlinear(out, args::Vararg; g::Function, g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(Nonlinear{Unscented}))
    return Nonlinear{Unscented}(out, args...; g=g, g_inv=g_inv, alpha=alpha, dims=dims, id=id)
end

slug(::Type{Nonlinear}) = "g"
