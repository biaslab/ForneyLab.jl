export Nonlinear, NonlinearUT, NonlinearPT, applyUnscentedTransform, applyParticleTransform

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

    function Nonlinear(out, in1, g::Function; id=ForneyLab.generateId(Nonlinear))
        @ensureVariables(out, in1)
        self = new(id, Vector{Interface}(undef, 2), Dict{Symbol,Interface}(), g)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
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
        graph = currentGraph()
        
        graph.nodes[node.id] = self
        node.i[:out].node = self
        node.i[:in1].node = self

        self.i[:out] = self.interfaces[1] = node.i[:out]
        self.i[:in1] = self.interfaces[2] = node.i[:in1]

        return self
    end
end

slug(::Type{NonlinearUT}) = "g"

function applyUnscentedTransform(edge::Edge; g_inv=nothing, alpha=nothing, dims=())
    applied = false
    if (edge.a !== nothing) && (edge.a.node isa Nonlinear)
        edge.a.node = NonlinearUT(edge.a.node, g_inv=g_inv, alpha=alpha, dims=dims)
        applied = true
    end
    if (edge.b !== nothing) && (edge.b.node isa Nonlinear)
        edge.b.node = NonlinearUT(edge.b.node, g_inv=g_inv, alpha=alpha, dims=dims)
        applied = true
    end
    return applied
end

function applyUnscentedTransform(var::Variable; g_inv=nothing, alpha=nothing, dims=())
    # find connected Nonlinear node(s)
    applied = false
    for edge in edges(var)
        applied = applied || applyUnscentedTransform(edge, g_inv=g_inv, alpha=alpha, dims=dims)
    end

    if !applied
        error("Cannot apply unscented transform to variable $(var.id).")
    end
end

function applyUnscentedTransform(vars::Vector{Variable}; g_inv=nothing, alpha=nothing, dims=())
    for var in vars
        applyUnscentedTransform(var, g_inv=g_inv, alpha=alpha, dims=dims)
    end
end

function applyUnscentedTransform(;g_inv=nothing, alpha=nothing, dims=())
    for node in currentGraph.nodes
        if node isa Nonlinear
            NonlinearUT(node, g_inv=g_inv, alpha=alpha, dims=dims)
        end
    end
end

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
        graph = currentGraph()
        
        graph.nodes[node.id] = self
        node.i[:out].node = self
        node.i[:in1].node = self

        self.i[:out] = self.interfaces[1] = node.i[:out]
        self.i[:in1] = self.interfaces[2] = node.i[:in1]
        
        return self
    end
end

slug(::Type{NonlinearPT}) = "g"

function applyParticleTransform(edge::Edge)
    applied = false
    if (edge.a !== nothing) && (edge.a.node isa Nonlinear)
        edge.a.node = NonlinearPT(edge.a.node)
        applied = true
    end
    if (edge.b !== nothing) && (edge.b.node isa Nonlinear)
        edge.b.node = NonlinearPT(edge.b.node)
        applied = true
    end
    return applied
end

function applyParticleTransform(var::Variable)
    # find connected Nonlinear node(s)
    applied = false
    for edge in edges(var)
        applied = applied || applyParticleTransform(edge)
    end

    if !applied
        error("Cannot apply particle transform to variable $(var.id).")
    end
end

function applyParticleTransform(vars::Vector{Variable})
    for var in vars
        unscentedTransform(var)
    end
end

function applyParticleTransform()
    for node in currentGraph.nodes
        if node isa Nonlinear
            node = NonlinearPT(node)
        end
    end
end
