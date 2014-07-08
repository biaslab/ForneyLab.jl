############################################
# AdditionNode
############################################
# Description:
#   Addition of two messages of same type.
#
#          in2
#          |
#    in1   v  out
#   ----->[+]----->
#
#   out = in1 + in2
#
#   Example:
#       AdditionNode(; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       GaussianMessage
#       GeneralMessage
#   2. (in2):
#       GaussianMessage
#       GeneralMessage
#   3. (out):
#       GaussianMessage
#       GeneralMessage
############################################

export AdditionNode

type AdditionNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    in2::Interface
    out::Interface

    function AdditionNode(; name="unnamed", args...)
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
# GaussianMessage methods
############################################

# Rule set for forward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardAdditionMRule{T<:Number}(m_x::Array{T, 1}, m_y::Array{T, 1}) = m_x + m_y
forwardAdditionVRule{T<:Number}(V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x + V_y
forwardAdditionWRule{T<:Number}(W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x * pinv(W_x + W_y) * W_y
forwardAdditionXiRule{T<:Number}(V_x::Array{T, 2}, xi_x::Array{T, 1}, V_y::Array{T, 2}, xi_y::Array{T, 1}) = pinv(V_x + V_y) * (V_x*xi_x + V_y*xi_y)

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
# The backward propagation merely negates the mean of the present input message (edge X) and uses the same rules to determine the missing input (edge Y)
# For the sake of clarity there is some redundancy between forward and backward rules.
backwardAdditionMRule{T<:Number}(m_x::Array{T, 1}, m_z::Array{T, 1}) = m_z - m_x
backwardAdditionVRule{T<:Number}(V_x::Array{T, 2}, V_z::Array{T, 2}) = V_x + V_z
backwardAdditionWRule{T<:Number}(W_x::Array{T, 2}, W_z::Array{T, 2}) = W_x * pinv(W_x + W_z) * W_z
backwardAdditionXiRule{T<:Number}(V_x::Array{T, 2}, xi_x::Array{T, 1}, V_z::Array{T, 2}, xi_z::Array{T, 1}) = pinv(V_x + V_z) * (V_z*xi_z - V_x*xi_x)

function updateNodeMessage!(outbound_interface_id::Int,
                            node::AdditionNode,
                            inbound_messages_types::Type{GaussianMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    msg_out = getOrAssign(node.interfaces[outbound_interface_id], GaussianMessage)

    # Calculations for the GaussianMessage type; Korl (2005), table 4.1
    if outbound_interface_id == 3
        # Forward message, both messages on the incoming edges, required to calculate the outgoing message.
        msg_1 = node.interfaces[1].partner.message
        msg_2 = node.interfaces[2].partner.message

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1.m != nothing && msg_1.V != nothing && msg_2.m != nothing && msg_2.V != nothing
            msg_out.m = forwardAdditionMRule(msg_1.m, msg_2.m)
            msg_out.V = forwardAdditionVRule(msg_1.V, msg_2.V)
            msg_out.W = nothing
            msg_out.xi= nothing
        elseif msg_1.m != nothing && msg_1.W != nothing && msg_2.m != nothing && msg_2.W != nothing
            msg_out.m = forwardAdditionMRule(msg_1.m, msg_2.m)
            msg_out.V = nothing
            msg_out.W = forwardAdditionWRule(msg_1.W, msg_2.W)
            msg_out.xi= nothing
        elseif msg_1.xi != nothing && msg_1.V != nothing && msg_2.xi != nothing && msg_2.V != nothing
            msg_out.m = nothing
            msg_out.V = forwardAdditionVRule(msg_1.V, msg_2.V)
            msg_out.W = nothing
            msg_out.xi= forwardAdditionXiRule(msg_1.V, msg_1.xi, msg_2.V, msg_2.xi)
        else
            # Last resort: calculate (m,V) parametrization for both inbound messages
            ensureMVParametrization!(msg_1)
            ensureMVParametrization!(msg_2)
            msg_out.m = forwardAdditionMRule(msg_1.m, msg_2.m)
            msg_out.V = forwardAdditionVRule(msg_1.V, msg_2.V)
            msg_out.W = nothing
            msg_out.xi= nothing
        end
    elseif outbound_interface_id == 1 || outbound_interface_id == 2
        # Backward message, one message on the incoming edge and one on the outgoing edge.
        msg_1or2 = (outbound_interface_id==1) ? node.interfaces[2].partner.message : node.interfaces[1].partner.message
        msg_3 = node.interfaces[3].partner.message

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1or2.m != nothing && msg_1or2.V != nothing && msg_3.m != nothing && msg_3.V != nothing
            msg_out.m = backwardAdditionMRule(msg_1or2.m, msg_3.m)
            msg_out.V = backwardAdditionVRule(msg_1or2.V, msg_3.V)
            msg_out.W = nothing
            msg_out.xi = nothing
        elseif msg_1or2.m != nothing && msg_1or2.W != nothing && msg_3.m != nothing && msg_3.W != nothing
            msg_out.m = backwardAdditionMRule(msg_1or2.m, msg_3.m)
            msg_out.V = nothing
            msg_out.W = backwardAdditionWRule(msg_1or2.W, msg_3.W)
            msg_out.xi = nothing
        elseif msg_1or2.xi != nothing && msg_1or2.V != nothing && msg_3.xi != nothing && msg_3.V != nothing
            msg_out.m = nothing
            msg_out.V = backwardAdditionVRule(msg_1or2.V, msg_3.V)
            msg_out.W = nothing
            msg_out.xi = backwardAdditionXiRule(msg_1or2.V, msg_1or2.xi, msg_3.V, msg_3.xi)
        else
            # Last resort: calculate (m,V) parametrization for both inbound messages
            ensureMVParametrization!(msg_1or2)
            ensureMVParametrization!(msg_3)
            msg_out.m = backwardAdditionMRule(msg_1or2.m, msg_3.m)
            msg_out.V = backwardAdditionVRule(msg_1or2.V, msg_3.V)
            msg_out.W = nothing
            msg_out.xi = nothing
        end
    else
        error("Invalid interface id ", outbound_interface_id, " for calculating message on ", typeof(node), " ", node.name)
    end

    return msg_out
end

#############################################
# GeneralMessage methods
#############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::AdditionNode,
                            inbound_messages_types::Type{GeneralMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    msg_out = getOrAssign(node.interfaces[outbound_interface_id], GeneralMessage)

    # Calculations for a general message type
    if outbound_interface_id == 1
        # Backward message 1
        msg_out.value = node.interfaces[3].partner.message.value - node.interfaces[2].partner.message.value
    elseif outbound_interface_id == 2
        # Backward message 2
        msg_out.value = node.interfaces[3].partner.message.value - node.interfaces[1].partner.message.value
    elseif outbound_interface_id == 3
        # Forward message
        msg_out.value = node.interfaces[1].partner.message.value + node.interfaces[2].partner.message.value
    else
        error("Invalid interface id ", outbound_interface_id, " for calculating message on ", typeof(node), " ", node.name)
    end

    return msg_out
end