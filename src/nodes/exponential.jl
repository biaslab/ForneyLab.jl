export ExponentialNode

"""
Description:

    Maps a Gaussian to a log-normal distribution.
    Derivations can be found in the derivations document.

    in         out
    ----->[exp]----->

    f(in,out) = δ(out - exp(in))

Interfaces:

    1 i[:in], 2 i[:out]

Construction:

    ExponentialNode(id=:my_node)
"""
type ExponentialNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function ExponentialNode(; id=generateNodeId(ExponentialNode))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:in, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::ExponentialNode) = true


############################################
# Exact Gaussian update functions
############################################

"""
ExponentialNode:

     N       logN
    --->[exp]--->
              -->
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::LogNormal,
                            msg_in::Message{Gaussian},
                            msg_out::Any)

    ensureParameters!(msg_in.payload, (:m, :V))

    outbound_dist.m = msg_in.payload.m
    outbound_dist.s = msg_in.payload.V

    return outbound_dist
end

"""
ExponentialNode:

     N       logN
    --->[exp]--->
    <--
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            msg_in::Any,
                            msg_out::Message{LogNormal})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.V = msg_out.payload.s
    outbound_dist.W = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end

############################################
# Approximate Gaussian update functions
############################################

"""
ExponentialNode:

     N       Gam (MomentMatching)
    --->[exp]--->
              -->

    The gamma msg is the moment-matching approx. to the exact log-normal message.
    The approximation is easily obtained by substituting the expressions for the mean and variance.
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Gamma,
                            msg_in::Message{Gaussian},
                            msg_out::Any,
                            approx::Type{MomentMatching})

    ensureParameters!(msg_in.payload, (:m, :V))
    outbound_dist.b = msg_in.payload.m / msg_in.payload.V
    outbound_dist.a = msg_in.payload.m * outbound_dist.b

    return outbound_dist
end

"""
ExponentialNode:

     N       Gam (LogMomentMatching)
    --->[exp]--->
              -->

    The approximate gamma msg matches the moments of ln(X).
    (The incoming Gaussian is approximated by p(ln X), where X is gamma distributed.)
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Gamma,
                            msg_in::Message{Gaussian},
                            msg_out::Any,
                            approx::Type{LogMomentMatching})

    ensureParameters!(msg_in.payload, (:m, :V))
    outbound_dist.a = trigammaInverse(msg_in.payload.V)
    outbound_dist.b = 1 / exp(msg_in.payload.m - digamma(outbound_dist.a))

    return outbound_dist
end

"""
ExponentialNode:

     N       Gam
    --->[exp]--->
    <--

    The Gaussian msg is the moment-matching approximation to the exact outbound message.
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            msg_in::Any,
                            msg_out::Message{Gamma},
                            approx::Type{MomentMatching})

    outbound_dist.m = digamma(msg_out.payload.a) - log(msg_out.payload.b)
    outbound_dist.V = trigamma(msg_out.payload.a)
    outbound_dist.W = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end

############################################
# Delta update functions
############################################

"""
ExponentialNode:

     δ        δ
    --->[exp]--->
              -->
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Delta{Float64},
                            msg_in::Message{Delta{Float64}},
                            msg_out::Any)

    outbound_dist.m = exp(msg_in.payload.m)
    return outbound_dist
end

"""
ExponentialNode:

     δ        δ
    --->[exp]--->
    <--
"""
function sumProductRule!(   node::ExponentialNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Delta{Float64},
                            msg_in::Any,
                            msg_out::Message{Delta{Float64}})

    outbound_dist.m = log(msg_out.payload.m)
    return outbound_dist
end


############################################
# MvGaussian update functions
############################################

function sumProductRule!{dims}( node::ExponentialNode,
                                outbound_interface_index::Type{Val{2}},
                                outbound_dist::MvLogNormal{dims},
                                msg_in::Message{MvGaussian{dims}},
                                msg_out::Any)

    ensureParameters!(msg_in.payload, (:m, :V))

    outbound_dist.m = deepcopy(msg_in.payload.m)
    outbound_dist.S = deepcopy(msg_in.payload.V)

    return outbound_dist
end

function sumProductRule!{dims}( node::ExponentialNode,
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussian{dims},
                                msg_in::Any,
                                msg_out::Message{MvLogNormal{dims}})

    outbound_dist.m = deepcopy(msg_out.payload.m)
    outbound_dist.V = deepcopy(msg_out.payload.S)
    invalidate!(outbound_dist.W)
    invalidate!(outbound_dist.xi)

    return outbound_dist
end


############################################
# MvDelta update functions
############################################

function sumProductRule!{T<:MvDelta{Float64}}(  node::ExponentialNode,
                                                            outbound_interface_index::Type{Val{2}},
                                                            outbound_dist::T,
                                                            msg_in::Message{T},
                                                            msg_out::Any)

    outbound_dist.m = exp(msg_in.payload.m)
    return outbound_dist
end

function sumProductRule!{T<:MvDelta{Float64}}(  node::ExponentialNode,
                                                            outbound_interface_index::Type{Val{1}},
                                                            outbound_dist::T,
                                                            msg_in::Any,
                                                            msg_out::Message{T})

    outbound_dist.m = log(msg_out.payload.m)
    return outbound_dist
end
