export Bivariate, Laplace

abstract type ApproximationMethod end
abstract type Laplace <: ApproximationMethod end

"""
Description:

    Bivariate node models a function of two random variables.

Interfaces:

    1. out
    2. in1
    3. in2

Construction:

    Bivariate(out, in1, in2, g, id=:my_node)
"""
mutable struct Bivariate{T<:ApproximationMethod} <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    dims::Tuple # Dimension of breaker message on input interface
    status::Dict{Symbol, Union{Integer, Message}} #Keeps the status of node to ensure that input variables are updated simultaneously

    function Bivariate{Laplace}(out, in1, in2, g::Function; id=ForneyLab.generateId(Bivariate{Laplace}))
        @ensureVariables(out, in1, in2)
        self = new(id, Vector{Interface}(undef, 3), Dict{Symbol,Interface}(), g, (), Dict(:count_update=>0))
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:in2] = self.interfaces[3] = associate!(Interface(self), in2)

        return self
    end
end

slug(::Type{Bivariate}) = "g"