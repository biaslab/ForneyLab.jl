export Nonlinear, Unscented, Sampling, Extended

abstract type ApproximationMethod end
abstract type Unscented <: ApproximationMethod end
abstract type Sampling <: ApproximationMethod end
abstract type Extended <: ApproximationMethod end

"""
Description:

    Nonlinear node modeling a nonlinear relation. Updates for
    the nonlinear node are computed through the unscented transform (by default), 
    importance sampling, or local linear approximation.

    For more details see "On Approximate Nonlinear Gaussian Message Passing on
    Factor Graphs", Petersen et al. 2018.

    f(out, in1) = δ(out - g(in1))

Interfaces:

    1. out
    2. in1

Construction:

    Nonlinear{T}(out, in1; g=g, id=:my_node)
    Nonlinear{T}(out, in1; g=g, g_inv=g_inv, id=:my_node)
    Nonlinear{T}(out, in1, in2, ...; g=g, id=:my_node)
    Nonlinear{T}(out, in1, in2, ...; g=g, g_inv=(g_inv_in1, g_inv_in2, ...), id=:my_node)
    Nonlinear{T}(out, in1, in2, ...; g=g, g_inv=(g_inv_in1, nothing, ...), id=:my_node)

    where T encodes the approximation method: Unscented, Sampling, or Extended.
"""
mutable struct Nonlinear{T<:ApproximationMethod} <: DeltaFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output as a function of the inputs
    g_inv::Union{Function, Nothing, Vector} # Inverse of g with respect to individual inbounds (optional)
    alpha::Union{Float64, Nothing} # Spread parameter for unscented transform
    dims::Union{Tuple, Vector} # Dimension of breaker message(s) on input interface(s)
    n_samples::Union{Int64, Nothing} # Number of samples for sampling

    function Nonlinear{T}(id::Symbol,
                          g::Function,
                          g_inv::Union{Function, Nothing, Vector},
                          alpha::Union{Float64, Nothing},
                          dims::Union{Tuple, Vector},
                          n_samples::Union{Int64, Nothing},
                          out::Variable,
                          args::Vararg) where T<:ApproximationMethod
        
        n_args = length(args)
        for i=1:n_args
            @ensureVariables(args[i])
        end
        self = new(id, Vector{Interface}(undef, n_args+1), Dict{Symbol,Interface}(), g, g_inv, alpha, dims, n_samples)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        for k = 1:n_args
            self.i[:in*k] = self.interfaces[k+1] = associate!(Interface(self), args[k])
        end

        return self
    end
end

function Nonlinear{Unscented}(out::Variable, args::Vararg; g::Function, g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(Nonlinear{Unscented}))
    return Nonlinear{Unscented}(id, g, g_inv, alpha, dims, nothing, out, args...)
end

function Nonlinear{Sampling}(out::Variable, args::Vararg; g::Function, dims=(), n_samples=nothing, id=ForneyLab.generateId(Nonlinear{Sampling}))
    return Nonlinear{Sampling}(id, g, nothing, nothing, dims, n_samples, out, args...)
end

function Nonlinear{Extended}(out::Variable, args::Vararg; g::Function, g_inv=nothing, dims=(), id=ForneyLab.generateId(Nonlinear{Extended}))
    return Nonlinear{Extended}(id, g, g_inv, nothing, dims, nothing, out, args...)
end

# A breaker message is required if interface is partnered with a Nonlinear{Sampling} inbound and there are multiple inbounds
function requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::Nonlinear{Sampling})
    backward = (partner_interface != partner_node.i[:out]) # Interface is partnered with an inbound
    multi_in = (length(partner_node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node

    return backward && multi_in
end

# A breaker message is required if interface is partnered with a Nonlinear{Unscented} inbound and no inverse is available
function requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::Nonlinear{T}) where T<:Union{Unscented, Extended}
    backward = (partner_interface != partner_node.i[:out]) # Interface is partnered with an inbound
    multi_in = (length(partner_node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node
    inx = findfirst(isequal(partner_interface), partner_node.interfaces) - 1 # Find number of inbound interface; 0 for outbound
    undefined_inverse = (partner_node.g_inv == nothing) || (multi_in && (inx > 0) && (partner_node.g_inv[inx] == nothing)) # (Inbound-specific) inverse is undefined

    return backward && undefined_inverse
end

function breakerParameters(interface::Interface, partner_interface::Interface, partner_node::Nonlinear)
    requiresBreaker(interface, partner_interface, partner_node) || error("Breaker dimensions requested for non-breaker interface: $(interface)")

    if isa(partner_node.dims, Vector) # Varying inbound dimensionalities are defined
        inx = findfirst(isequal(partner_interface), partner_node.interfaces) - 1 # Find number of inbound interface
        dims = partner_node.dims[inx] # Extract dimensionality from vector
    else
        dims = partner_node.dims # All inbound dimensions are equal
    end

    # Deterimine message variate type
    if dims == ()
        var_type = Univariate
    else
        var_type = Multivariate
    end

    return (Message{GaussianMeanVariance, var_type}, dims)
end

slug(::Type{Nonlinear{T}}) where T<:ApproximationMethod = "g{$(removePrefix(T))}"
