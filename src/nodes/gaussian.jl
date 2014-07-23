############################################
# GaussianNode
############################################
# Description:
#   Node converting an input mean and precision
#   to a univariate Gaussian distribution.
#   
#   The GaussianNode has two different names for its second interface,
#   namely precision and variance. These named handles are simply
#   pointers to one and the same interface. 
#
#         mean
#          |
#          v  out
#   ----->[N]----->
#  precision/
#  variance
#
#   out = Message(GaussianDistribution(m=mean, W=precision))
#
#   Example:
#       GaussianNode([1.0], [0.1]; name="my_node")
#
# Interface ids, (names) and supported message types:
#   Receiving:
#   1. (mean):
#       Message{Float64}
#       GaussianDistribution (marginal)
#   2. (precision / variance):
#       Message{Float64}
#       GammaDistribution (marginal)
#       InverseGammaDistribution (marginal)
#   3. (out):
#       Message{Float64}
#
#   Sending:
#   1. (mean):
#       Message{GaussianDistribution}
#   2. (precision / variance):
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   3. (out):
#       Message{GaussianDistribution}
############################################

export GaussianNode

type GaussianNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    variational::Bool
    # Helper fields filled by constructor
    mean::Interface
    precision::Interface # alias
    variance::Interface # alias
    out::Interface
    function GaussianNode(variational; name="unnamed", args...)
        self = new(name, Array(Interface, 3), variational)
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Init named interface handles
        self.mean = self.interfaces[1]
        self.precision = self.interfaces[2] # alias
        self.variance = self.interfaces[2] # alias
        self.out = self.interfaces[3]
        return self
    end
end
GaussianNode(; args...) = GaussianNode(false; args...)

############################################
# Update functions
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::GaussianNode,
                            inbound_messages_value_types::Any,
                            outbound_message_value_type::Any)
    # Select update function depending on the variational or observed nature of the Gaussian node
    if node.variational
        updateNodeMessageVariational!(outbound_interface_id, node, inbound_messages_value_types, outbound_message_value_type)
    else
        updateNodeMessagePointEstimate!(outbound_interface_id, node, inbound_messages_value_types, outbound_message_value_type)
    end
end

############################################
# Point-estimate update functions
############################################

function updateNodeMessagePointEstimate!(outbound_interface_id::Int,
                                         node::GaussianNode,
                                         inbound_messages_value_types::Type{Union(Float64, InverseGammaDistribution)},
                                         outbound_message_value_type::Type{GaussianDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    # For interface 1 the calling signature for this function is identical to the variational call.
    # A decision about where to direct the call is made by the invoking updateNodeMessage function.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Formulas from table 5.2 in Korl (2005)
    if outbound_interface_id == 3
        # Forward message
        gamma = node.interfaces[2].partner.message.value
        mean = node.interfaces[1].partner.message.value
        dist_out.m = [mean]
        dist_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    elseif outbound_interface_id == 1
        # Backward over mean edge
        # Rules not in Korl, but equivalent by symmetry
        gamma = node.interfaces[2].partner.message.value
        y = node.interfaces[3].partner.message.value
        dist_out.m = [y]
        dist_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end
function updateNodeMessagePointEstimate!(outbound_interface_id::Int,
                                         node::GaussianNode,
                                         inbound_messages_value_types::Type{Union(Float64, GammaDistribution)},
                                         outbound_message_value_type::Type{GaussianDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    error("Point estimate update rules for gamma message on Gaussian node not implemented yet")
end

function updateNodeMessagePointEstimate!(outbound_interface_id::Int,
                                         node::GaussianNode,
                                         inbound_messages_value_types::Type{Float64},
                                         outbound_message_value_type::Type{InverseGammaDistribution})

    # Calculate the message for the variance when both m and y are incoming.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Formulas from table 5.2 in Korl (2005)
    if outbound_interface_id == 2
        # Backward over variance edge
        y = node.interfaces[3].partner.message.value
        m = node.interfaces[1].partner.message.value
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end
    # Set the outbound message
    return node.interfaces[outbound_interface_id].message
end

############################################
# Variational update functions
############################################

function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_value_types::Type{Union(Float64, InverseGammaDistribution)},
                                       outbound_message_value_type::Type{GaussianDistribution})

    # Variational update function, takes the marginals as input (instead of the inbound messages)
    # and update the outgoing message on the mean interface of a Gaussian node.
    # For interface 1 the calling signature for this function is identical to the standard sumproduct call.
    # A decision about where to direct the call is made by the invoking updateNodeMessage function.
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if outbound_interface_id == 1 # Mean estimation from variance and sample
        y_0 = node.out.partner.message.value # observation
        a = node.variance.edge.marginal.a # gamma message
        b = node.variance.edge.marginal.b
        dist_out.m = [y_0]
        dist_out.V = reshape([((a+1)/b)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return node.interfaces[outbound_interface_id].message
end
function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_value_types::Type{Union(Float64, GammaDistribution)},
                                       outbound_message_value_type::Type{GaussianDistribution})

    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the mean interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if outbound_interface_id == 1 # Mean estimation from variance and sample
        y_0 = node.out.partner.message.value # observation
        a = node.precision.edge.marginal.a # gamma distribution
        b = node.precision.edge.marginal.b
        dist_out.m = [y_0]
        dist_out.V = nothing
        dist_out.xi = nothing
        dist_out.W = reshape([a/b], 1, 1)
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_value_types::Type{Union(Float64, GaussianDistribution)},
                                       outbound_message_value_type::Type{GammaDistribution})
    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the variance interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    if outbound_interface_id == 2 # Variance/precision estimation from mean and sample
        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value
        y_0 = node.out.partner.message.value # observation
        ensureMWParametrization!(node.mean.edge.marginal)
        m = node.mean.edge.marginal.m[1] # Gaussian distribution
        W = node.mean.edge.marginal.W[1,1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y_0-m)^2+0.5*inv(W)
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_value_types::Type{Union(Float64, GaussianDistribution)},
                                       outbound_message_value_type::Type{InverseGammaDistribution})
    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the variance interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    if outbound_interface_id == 2 # Variance/precision estimation from mean and sample
        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value
        y_0 = node.out.partner.message.value # observation
        ensureMVParametrization!(node.mean.edge.marginal)
        m = node.mean.edge.marginal.m[1] # Gaussian message
        V = node.mean.edge.marginal.V[1,1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y_0-m)^2+0.5*V
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end