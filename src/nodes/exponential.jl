############################################
# ExponentialNode
############################################
# Description:
#   Maps a Gaussian to a log-normal distribution.
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
function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Int,
                        msg_in::Message{GaussianDistribution},
                        msg_out::Void)
    dist_out = ensureMessage!(node.i[:out], LogNormalDistribution).payload

    ensureParameters!(msg_in.payload, (:m, :V))

    dist_out.m = msg_in.payload.m
    dist_out.s = msg_in.payload.V

    return (:exponential_forward_gaussian,
            node.interfaces[outbound_interface_index].message)
end

# Backward message
function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Int,
                        msg_in::Void,
                        msg_out::Message{LogNormalDistribution})
    dist_out = ensureMessage!(node.i[:in], GaussianDistribution).payload

    dist_out.m = msg_out.payload.m
    dist_out.V = msg_out.payload.s
    dist_out.W = NaN
    dist_out.xi = NaN

    return (:exponential_backward_gaussian,
            node.interfaces[outbound_interface_index].message)
end


############################################
# DeltaDistribution update functions
############################################

# Forward message
function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Int,
                        msg_in::Message{DeltaDistribution{Float64}},
                        msg_out::Void)
    dist_out = ensureMessage!(node.i[:out], DeltaDistribution{Float64}).payload

    dist_out.m = exp(msg_in.payload.m)

    return (:exponential_forward_delta,
            node.interfaces[outbound_interface_index].message)
end

# Backward message
function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Int,
                        msg_in::Void,
                        msg_out::Message{DeltaDistribution{Float64}})
    dist_out = ensureMessage!(node.i[:in], DeltaDistribution{Float64}).payload

    dist_out.m = log(msg_out.payload.m)

    return (:exponential_backward_delta,
            node.interfaces[outbound_interface_index].message)
end
