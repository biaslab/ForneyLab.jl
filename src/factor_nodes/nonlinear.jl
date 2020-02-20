export Nonlinear, Unscented, ImportanceSampling

abstract type ApproximationMethod end
abstract type Unscented <: ApproximationMethod end
abstract type ImportanceSampling <: ApproximationMethod end

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

    Nonlinear(out, in1, g, id=:my_node)
"""
mutable struct Nonlinear{T<:ApproximationMethod} <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    g_inv::Union{Function, Nothing} # Inverse of g (optional)
    alpha::Union{Float64, Nothing} # Spread parameter for unscented transform
    dims::Tuple # Dimension of breaker message on input interface

    function Nonlinear{Unscented}(out, in1, g::Function; g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(Nonlinear{Unscented}))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), g, g_inv, alpha, dims)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end

    function Nonlinear{ImportanceSampling}(out, in1, g::Function; id=ForneyLab.generateId(Nonlinear{ImportanceSampling}))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), g, nothing, nothing, ())
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

function Nonlinear(out, in1, g::Function; g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(Nonlinear{Unscented}))
    return Nonlinear{Unscented}(out, in1, g::Function; g_inv=g_inv, alpha=alpha, dims=dims, id=id)
end

slug(::Type{Nonlinear}) = "g"

