############################################
# LogarithmicNode
############################################
# Description:
#   Maps a Gaussian to a gamma distribution.
#   Derivations can be found in the derivations document.
#
#    in1        out
#   ----->[log]----->
#
#   Example:
#       LogarithmicNode(; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       Message{GaussianDistribution}
#       Message{DeltaDistribution}
#   2. (out):
#       Message{GammaDistribution}
#       Message{DeltaDistribution}
############################################

export LogarithmicNode

type LogarithmicNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    out::Interface

    function LogarithmicNode(; name=unnamedStr())
        self = new(name, Array(Interface, 2))

        named_handle_list = [:in1, :out]
        for i = 1:length(named_handle_list)
            self.interfaces[i] = Interface(self)
            setfield!(self, named_handle_list[i], self.interfaces[i])
        end

        return self
    end
end

isDeterministic(::LogarithmicNode) = true


############################################
# Standard update functions
############################################

# Forward message
function sumProduct!(node::LogarithmicNode,
                     outbound_interface_id::Int,
                     msg_in1::Message{GaussianDistribution},
                     msg_out::Nothing)
    dist_out = getOrCreateMessage(node.out, GammaDistribution).payload

    ensureMWParametrization!(msg_in1.payload)
    (length(msg_in1.payload.m) == 1) || error("Forward update for LogarithmicNode only defined for univariate input")

    gam = msg_in1.payload.W[1,1]
    mu = msg_in1.payload.m[1]

    dist_out.a = gam + 1
    dist_out.b = gam/(exp(mu))

    return node.interfaces[outbound_interface_id].message
end

# Backward message
function sumProduct!(node::LogarithmicNode,
                     outbound_interface_id::Int,
                     msg_in1::Nothing,
                     msg_out::Message{GammaDistribution})
    dist_out = getOrCreateMessage(node.in1, GaussianDistribution).payload

    a = msg_out.payload.a
    b = msg_out.payload.b

    dist_out.m = [log((a-1)/b)]
    dist_out.V = nothing
    dist_out.W = reshape([a-1], 1, 1)
    dist_out.xi= nothing

    return node.interfaces[outbound_interface_id].message
end


############################################
# DeltaDistribution update functions
############################################

# Forward message
function sumProduct!{T<:Any}(node::LogarithmicNode,
                     outbound_interface_id::Int,
                     msg_in1::Message{DeltaDistribution{T}},
                     msg_out::Nothing)
    length(msg_in1.payload.m) == 1 || error("LogarithmicNode only defined for univariate variables")
    dist_out = getOrCreateMessage(node.out, DeltaDistribution{T}).payload

    dist_out.m = log(msg_in1.payload.m)

    return node.interfaces[outbound_interface_id].message
end

# Backward message
function sumProduct!{T<:Any}(node::LogarithmicNode,
                     outbound_interface_id::Int,
                     msg_in1::Nothing,
                     msg_out::Message{DeltaDistribution{T}})
    length(msg_out.payload.m) == 1 || error("LogarithmicNode only defined for univariate variables")
    dist_out = getOrCreateMessage(node.in1, DeltaDistribution{T}).payload

    dist_out.m = exp(msg_out.payload.m)

    return node.interfaces[outbound_interface_id].message
end
