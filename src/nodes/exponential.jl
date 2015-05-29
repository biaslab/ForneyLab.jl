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
#   ExponentialNode(name="my_node")
#
############################################

export ExponentialNode

type ExponentialNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function ExponentialNode(; name=unnamedStr())
        self = new(name, Array(Interface, 2))

        for (iface_id, iface_name) in enumerate([:in, :out])
            self.i[iface_name] = self.interfaces[iface_id] = Interface(self)
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
                     outbound_interface_id::Int,
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
            node.interfaces[outbound_interface_id].message)
end

# Backward message
function sumProduct!(node::ExponentialNode,
                     outbound_interface_id::Int,
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
            node.interfaces[outbound_interface_id].message)
end


############################################
# DeltaDistribution update functions
############################################

# Forward message
function sumProduct!{T<:Any}(node::ExponentialNode,
                     outbound_interface_id::Int,
                     msg_in::Message{DeltaDistribution{T}},
                     msg_out::Nothing)
    length(msg_in.payload.m) == 1 || error("ExponentialNode only defined for univariate variables")
    dist_out = ensureMessage!(node.i[:out], DeltaDistribution{T}).payload

    dist_out.m = exp(msg_in.payload.m)

    return (:exponential_forward_delta,
            node.interfaces[outbound_interface_id].message)
end

# Backward message
function sumProduct!{T<:Any}(node::ExponentialNode,
                     outbound_interface_id::Int,
                     msg_in::Nothing,
                     msg_out::Message{DeltaDistribution{T}})
    length(msg_out.payload.m) == 1 || error("ExponentialNode only defined for univariate variables")
    dist_out = ensureMessage!(node.i[:in], DeltaDistribution{T}).payload

    dist_out.m = log(msg_out.payload.m)

    return (:exponential_backward_delta,
            node.interfaces[outbound_interface_id].message)
end
