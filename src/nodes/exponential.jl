############################################
# ExponentialNode
############################################
# Description:
#   Maps a Gaussian to a gamma distribution.
#   Derivations can be found in the derivations document.
#
#    in1        out
#   ----->[exp]----->
#
#   Example:
#       ExponentialNode(; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       Message{GaussianDistribution}
#       Message{DeltaDistribution}
#   2. (out):
#       Message{GammaDistribution}
#       Message{DeltaDistribution}
############################################

export ExponentialNode

type ExponentialNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    out::Interface

    function ExponentialNode(; name=unnamedStr())
        self = new(name, Array(Interface, 2))

        named_handle_list = [:in1, :out]
        for i = 1:length(named_handle_list)
            self.interfaces[i] = Interface(self)
            setfield!(self, named_handle_list[i], self.interfaces[i])
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
                     msg_in1::Message{GaussianDistribution},
                     msg_out::Nothing)
    dist_out = ensureMessage!(node.out, GammaDistribution).payload

    ensureMWParametrization!(msg_in1.payload)
    (length(msg_in1.payload.m) == 1) || error("Forward update for ExponentialNode only defined for univariate input")

    gam = msg_in1.payload.W[1,1]
    mu = msg_in1.payload.m[1]

    dist_out.a = gam + 1
    dist_out.b = gam/(exp(mu))

    return (:exponential_forward_gaussian,
            node.interfaces[outbound_interface_id].message)
end

# Backward message
function sumProduct!(node::ExponentialNode,
                     outbound_interface_id::Int,
                     msg_in1::Nothing,
                     msg_out::Message{GammaDistribution})
    dist_out = ensureMessage!(node.in1, GaussianDistribution).payload

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
                     msg_in1::Message{DeltaDistribution{T}},
                     msg_out::Nothing)
    length(msg_in1.payload.m) == 1 || error("ExponentialNode only defined for univariate variables")
    dist_out = ensureMessage!(node.out, DeltaDistribution{T}).payload

    dist_out.m = exp(msg_in1.payload.m)

    return (:exponential_forward_delta,
            node.interfaces[outbound_interface_id].message)
end

# Backward message
function sumProduct!{T<:Any}(node::ExponentialNode,
                     outbound_interface_id::Int,
                     msg_in1::Nothing,
                     msg_out::Message{DeltaDistribution{T}})
    length(msg_out.payload.m) == 1 || error("ExponentialNode only defined for univariate variables")
    dist_out = ensureMessage!(node.in1, DeltaDistribution{T}).payload

    dist_out.m = log(msg_out.payload.m)

    return (:exponential_backward_delta,
            node.interfaces[outbound_interface_id].message)
end
