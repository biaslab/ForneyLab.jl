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
#       GeneralMessage
#   2. (in2):
#       (Standard or inverted) GammaMessage
#   3. (out):
#       GeneralMessage
#
#   Sending:
#   1. (in1):
#       GaussianMessage
#   2. (in2):
#       (Inverted) GammaMessage
#   3. (out):
#       GaussianMessage
############################################

export GaussianNode

type GaussianNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface
    function GaussianNode(; args...)
        (name = getArgumentValue(args, :name))!=false || (name = "unnamed")
        self = new(name, Array(Interface, 3))
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

############################################
# Point-estimate update functions
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::GaussianNode,
                            inbound_messages_types::Type{Union(GeneralMessage, GammaMessage)})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    # Formulas from table 5.2 in Korl (2005)
    if outbound_interface_id == 3
        # Forward message
        gamma = node.interfaces[2].partner.message
        mean = node.interfaces[1].partner.message.value
        if gamma.inverted
            msg_out = GaussianMessage(m = [mean], V=[gamma.b/(gamma.a-1)]) # V is just the mode of the inverse gamma
        else
            msg_out = GaussianMessage(m = [mean], W=[(gamma.a-1)/gamma.b]) # W is just the mode of the gamma
        end
    elseif outbound_interface_id == 1

        # TODO: check duplicate function

        # Backward over mean edge
        # Rules not in Korl, but equivalent by symmetry
        gamma = node.interfaces[2].partner.message
        y = node.interfaces[3].partner.message.value
        if gamma.inverted
            msg_out = GaussianMessage(m = [y], V=[gamma.b/(gamma.a-1)]) # V is just the mode of the inverse gamma
        else
            msg_out = GaussianMessage(m = [y], W=[(gamma.a-1)/gamma.b]) # W is just the mode of the gamma
        end
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end
    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg_out
end
function updateNodeMessage!(outbound_interface_id::Int,
                            node::GaussianNode,
                            inbound_messages_types::Type{GeneralMessage})

    # Both m and y are known
    # Formulas from table 5.2 in Korl (2005)
    if outbound_interface_id == 2
        # Backward over variance edge
        y = node.interfaces[3].partner.message.value
        m = node.interfaces[1].partner.message.value
        msg_out = GammaMessage(a=-0.5, b=0.5*(y-m)^2, inverted=true) # Send inverse gamma message
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end
    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg_out
end

############################################
# Variational update functions
############################################

function updateNodeMessage!(outbound_interface_id::Int, node::GaussianNode, inbound_messages_types::Type{Union(GeneralMessage, GammaMessage)})
    
    # TODO: check duplicate function

    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the mean interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    if outbound_interface_id == 1 # Mean estimation from variance and sample
        y_0 = node.out.partner.message.value # observation
        a = node.in2.edge.marginal.a # gamma message
        b = node.in2.edge.marginal.b
        nu_m = GaussianMessage( m=[y_0], V=[((a+1)/b)] )
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return node.interfaces[outbound_interface_id].message = deepcopy(nu_m)
end
function updateNodeMessage!(outbound_interface_id::Int, node::GaussianNode, inbound_messages_types::Type{Union(GeneralMessage, GaussianMessage)})
    # Variational update function, takes the MARGINALS as input instead of the inbound messages.
    # Update the outgoing message on the variance interface of a Gaussian node.
    # Derivation for the update rule can be found in the derivations notebook.

    if outbound_interface_id == 2 # Variance estimation from mean and sample
        y_0 = node.out.partner.message.value # observation
        m = node.in1.edge.marginal.m[1] # Gaussian message
        V = node.in1.edge.marginal.V[1,1]
        nu_s = GammaMessage( a=-0.5, b=0.5*(y_0-m)^2+0.5*V, inverted=true )
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end
    return node.interfaces[outbound_interface_id].message = deepcopy(nu_s)
end
