############################################
# GaussianNode
############################################
# Description:
#   Node converting an input mean and precision
#   to a univariate Gaussian distribution.
#   Estimate precision with given inputs and vice versa,
#   Two distributions as input are not implemented 
#   and is a task for vmp:
#
#         in1 (mean)
#          |
#          v  out
#   ----->[N]----->
#  in2 (prec.)
#
#   out = GaussianMessage(m=in1, W=in2)
#
#   Example:
#       GaussianNode([1.0], [0.1]; name="my_node")
#
# Interface ids, (names) and supported message types:
#   Receiving:
#   1. (in1):
#       GeneralMessage / GaussianMessage marginal
#   2. (in2):
#       GeneralMessage / (Inverse)GammaMessage marginal
#   3. (out):
#       GeneralMessage
#
#   Sending:
#   1. (in1):
#       GaussianMessage
#   2. (in2):
#       (Inverse)GammaMessage
#   3. (out):
#       GaussianMessage
############################################

export GaussianNode

type GaussianNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    variational::Bool
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface
    function GaussianNode(variational; name="unnamed", args...)
        self = new(name, Array(Interface, 3), variational)
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.in2 = self.interfaces[2]
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
                            inbound_messages_types::Any)
    # Select update function depending on the variational or observed nature of the Gaussian node
    if node.variational
        updateNodeMessageVariational!(outbound_interface_id, node, inbound_messages_types)
    else
        updateNodeMessagePointEstimate!(outbound_interface_id, node, inbound_messages_types)
    end
end

############################################
# Point-estimate update functions
############################################

function updateNodeMessagePointEstimate!(outbound_interface_id::Int,
                                         node::GaussianNode,
                                         inbound_messages_types::Type{Union(GeneralMessage, InverseGammaMessage)})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    # For interface 1 this call is identical to the variational call. A decision about the right call is made by the invoking updateNodeMessage function.

    msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianMessage)

    # Formulas from table 5.2 in Korl (2005)
    if outbound_interface_id == 3
        # Forward message
        gamma = node.interfaces[2].partner.message
        mean = node.interfaces[1].partner.message.value
        msg_out.m = [mean]
        msg_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        msg_out.xi = nothing
        msg_out.W = nothing
    elseif outbound_interface_id == 1
        # Backward over mean edge
        # Rules not in Korl, but equivalent by symmetry
        gamma = node.interfaces[2].partner.message
        y = node.interfaces[3].partner.message.value
        msg_out.m = [y]
        msg_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        msg_out.xi = nothing
        msg_out.W = nothing
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end

    return msg_out
end
function updateNodeMessagePointEstimate!(outbound_interface_id::Int,
                                         node::GaussianNode,
                                         inbound_messages_types::Type{Union(GeneralMessage, GammaMessage)})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    # For interface 1 this call is identical to the variational call. A decision about the right call is made by the invoking updateNodeMessage function.

    error("Point estimate update rules for gamma message on Gaussian node not implemented yet")
end

function updateNodeMessagePointEstimate!(outbound_interface_id::Int,
                                         node::GaussianNode,
                                         inbound_messages_types::Type{GeneralMessage})

    # Both m and y are known

    msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], InverseGammaMessage)

    # Formulas from table 5.2 in Korl (2005)
    if outbound_interface_id == 2
        # Backward over variance edge
        y = node.interfaces[3].partner.message.value
        m = node.interfaces[1].partner.message.value
        msg_out.a = -0.5
        msg_out.b = 0.5*(y-m)^2
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end
    # Set the outbound message
    return msg_out
end

############################################
# Variational update functions
############################################

function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_types::Type{Union(GeneralMessage, InverseGammaMessage)})

    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # For interface 1 this call is identical to the point estimate call. A decision about the right call is made by the invoking updateNodeMessage function.
    # Update the outgoing message on the mean interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianMessage)

    if outbound_interface_id == 1 # Mean estimation from variance and sample
        y_0 = node.out.partner.message.value # observation
        a = node.in2.edge.marginal.a # gamma message
        b = node.in2.edge.marginal.b
        msg_out.m = [y_0]
        msg_out.V = reshape([((a+1)/b)], 1, 1)
        msg_out.xi = nothing
        msg_out.W = nothing
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return msg_out
end
function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_types::Type{Union(GeneralMessage, GammaMessage)})

    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the mean interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianMessage)

    if outbound_interface_id == 1 # Mean estimation from variance and sample
        y_0 = node.out.partner.message.value # observation
        a = node.in2.edge.marginal.a # gamma message
        b = node.in2.edge.marginal.b
        msg_out.m = [y_0]
        msg_out.V = nothing
        msg_out.xi = nothing
        msg_out.W = reshape([a/b], 1, 1)
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return msg_out
end

function updateNodeMessageVariational!(outbound_interface_id::Int,
                                       node::GaussianNode,
                                       inbound_messages_types::Type{Union(GeneralMessage, GaussianMessage)})
    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the variance interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    if outbound_interface_id == 2 # Variance/precision estimation from mean and sample
        if typeof(node.in2.edge.marginal) == GammaMessage # return standard gamma
            msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GammaMessage)
            y_0 = node.out.partner.message.value # observation
            ensureMWParametrization!(node.in1.edge.marginal)
            m = node.in1.edge.marginal.m[1] # Gaussian messsage
            W = node.in1.edge.marginal.W[1,1]
            msg_out.a = 1.5
            msg_out.b = 0.5*(y_0-m)^2+0.5*inv(W)
        elseif typeof(node.in2.edge.marginal) == InverseGammaMessage # Return inverse gamma 
            msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], InverseGammaMessage)
            y_0 = node.out.partner.message.value # observation
            ensureMVParametrization!(node.in1.edge.marginal)
            m = node.in1.edge.marginal.m[1] # Gaussian message
            V = node.in1.edge.marginal.V[1,1]
            msg_out.a = -0.5
            msg_out.b = 0.5*(y_0-m)^2+0.5*V
        else
            error("Don't know whether to return gamma or inverse gamma message. Please preset the marginal on the in2 interface's edge of GaussianNode $(node.name).")
        end
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return msg_out
end
