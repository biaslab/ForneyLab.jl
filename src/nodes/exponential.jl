############################################
# ExponentialNode
############################################
# Description:
#   Maps a Gaussian to a gamma distribution.
#   Derivations can be found in the derivations document.
#
#    in         out
#   ----->[exp]----->
#
#   f(in,out) = δ(out - exp(in))
#
# Interfaces:
#   1 i[:in], 2 i[:out]
#
# Construction:
#   ExponentialNode(id=:my_node)
#
############################################

export ExponentialNode

type ExponentialNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function ExponentialNode(; id=generateNodeId(ExponentialNode))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(current_graph, self)
 
        for (iface_index, iface_handle) in enumerate([:in, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::ExponentialNode) = true


############################################
# Standard update functions
############################################

# Forward message
function sumProduct!(node::ExponentialNode,
                     outbound_interface_index::Int,
                     msg_in::Message{GaussianDistribution},
                     msg_out::Nothing)
    dist_out = ensureMessage!(node.i[:out], GammaDistribution).payload

    ensureMWParametrization!(msg_in.payload)
    (length(msg_in.payload.m) == 1) || error("Forward update for ExponentialNode only defined for univariate input")

    gam = msg_in.payload.W[1,1]
    mu = msg_in.payload.m[1]

    dist_out.a = gam + 1
    dist_out.b = gam/(exp(mu))

    return (:exponential_forward_gaussian,
            node.interfaces[outbound_interface_index].message)
end

# Backward message
function sumProduct!(node::ExponentialNode,
                     outbound_interface_index::Int,
                     msg_in::Nothing,
                     msg_out::Message{GammaDistribution})
    dist_out = ensureMessage!(node.i[:in], GaussianDistribution).payload

    a = msg_out.payload.a
    b = msg_out.payload.b

    dist_out.m = [log((a-1)/b)]
    invalidate!(dist_out.V)
    dist_out.W = reshape([a-1], 1, 1)
    invalidate!(dist_out.xi)

    return (:exponential_backward_gaussian,
            node.interfaces[outbound_interface_index].message)
end


############################################
# DeltaDistribution update functions
############################################

# Forward message
function sumProduct!(node::ExponentialNode,
                     outbound_interface_index::Int,
                     msg_in::Message{DeltaDistribution{Float64}},
                     msg_out::Nothing)
    length(msg_in.payload.m) == 1 || error("ExponentialNode only defined for univariate variables")
    dist_out = ensureMessage!(node.i[:out], DeltaDistribution{Float64}).payload

    dist_out.m = exp(msg_in.payload.m)

    return (:exponential_forward_delta,
            node.interfaces[outbound_interface_index].message)
end

# Backward message
function sumProduct!(node::ExponentialNode,
                     outbound_interface_index::Int,
                     msg_in::Nothing,
                     msg_out::Message{DeltaDistribution{Float64}})
    length(msg_out.payload.m) == 1 || error("ExponentialNode only defined for univariate variables")
    dist_out = ensureMessage!(node.i[:in], DeltaDistribution{Float64}).payload

    dist_out.m = log(msg_out.payload.m)

    return (:exponential_backward_delta,
            node.interfaces[outbound_interface_index].message)
end
