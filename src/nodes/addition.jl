############################################
# AdditionNode
############################################
# Description:
#   Addition of two messages of same type.
#   out = in1 + in2
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

    function AdditionNode(;args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
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
                            inbound_messages::Array{GaussianMessage, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end

    # Calculations for the GaussianMessage type; Korl (2005), table 4.1
    if outbound_interface_id == 3
        # Forward message, both messages on the incoming edges, required to calculate the outgoing message.
        msg_in1 = inbound_messages[1]
        msg_in2 = inbound_messages[2]
        msg = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_in1.m != nothing && msg_in1.V != nothing && msg_in2.m != nothing && msg_in2.V != nothing
            msg.m = forwardAdditionMRule(msg_in1.m, msg_in2.m)
            msg.V = forwardAdditionVRule(msg_in1.V, msg_in2.V)
            msg.W = nothing
            msg.xi = nothing
        elseif msg_in1.m != nothing && msg_in1.W != nothing && msg_in2.m != nothing && msg_in2.W != nothing
            msg.m = forwardAdditionMRule(msg_in1.m, msg_in2.m)
            msg.V = nothing
            msg.W = forwardAdditionWRule(msg_in1.W, msg_in2.W)
            msg.xi = nothing
        elseif msg_in1.xi != nothing && msg_in1.V != nothing && msg_in2.xi != nothing && msg_in2.V != nothing
            msg.m = nothing
            msg.V = forwardAdditionVRule(msg_in1.V, msg_in2.V)
            msg.W = nothing
            msg.xi = forwardAdditionXiRule(msg_in1.V, msg_in1.xi, msg_in2.V, msg_in2.xi)
        else
            # Last resort: calculate (m,V) parametrization for both inbound messages
            ensureMVParametrization!(msg_in1)
            ensureMVParametrization!(msg_in2)
            msg.m = forwardAdditionMRule(msg_in1.m, msg_in2.m)
            msg.V = forwardAdditionVRule(msg_in1.V, msg_in2.V)
            msg.W = nothing
            msg.xi = nothing
        end
    elseif outbound_interface_id == 1 || outbound_interface_id == 2
        # Backward message, one message on the incoming edge and one on the outgoing edge.
        msg_in = (outbound_interface_id==1) ? inbound_messages[2] : inbound_messages[1]
        msg_out = inbound_messages[3]
        msg = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_in.m != nothing && msg_in.V != nothing && msg_out.m != nothing && msg_out.V != nothing
            msg.m = backwardAdditionMRule(msg_in.m, msg_out.m)
            msg.V = backwardAdditionVRule(msg_in.V, msg_out.V)
            msg.W = nothing
            msg.xi = nothing
        elseif msg_in.m != nothing && msg_in.W != nothing && msg_out.m != nothing && msg_out.W != nothing
            msg.m = backwardAdditionMRule(msg_in.m, msg_out.m)
            msg.V = nothing
            msg.W = backwardAdditionWRule(msg_in.W, msg_out.W)
            msg.xi = nothing
        elseif msg_in.xi != nothing && msg_in.V != nothing && msg_out.xi != nothing && msg_out.V != nothing
            msg.m = nothing
            msg.V = backwardAdditionVRule(msg_in.V, msg_out.V)
            msg.W = nothing
            msg.xi = backwardAdditionXiRule(msg_in.V, msg_in.xi, msg_out.V, msg_out.xi)
        else
            # Last resort: calculate (m,V) parametrization for both inbound messages
            ensureMVParametrization!(msg_in)
            ensureMVParametrization!(msg_out)
            msg.m = backwardAdditionMRule(msg_in.m, msg_out.m)
            msg.V = backwardAdditionVRule(msg_in.V, msg_out.V)
            msg.W = nothing
            msg.xi = nothing
        end
    else
        error("Invalid interface id ", outbound_interface_id, " for calculating message on ", typeof(node), " ", node.name)
    end

    return node.interfaces[outbound_interface_id].message = msg
end

#############################################
# GeneralMessage methods
#############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::AdditionNode,
                            inbound_messages::Array{GeneralMessage, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end

    # Calculations for a general message type
    if outbound_interface_id == 1
        # Backward message 1
        msg = GeneralMessage(inbound_messages[3].value - inbound_messages[2].value)
    elseif outbound_interface_id == 2
        # Backward message 2
        msg = GeneralMessage(inbound_messages[3].value - inbound_messages[1].value)
    elseif outbound_interface_id == 3
        # Forward message
        msg = GeneralMessage(inbound_messages[1].value + inbound_messages[2].value)
    else
        error("Invalid interface id ", outbound_interface_id, " for calculating message on ", typeof(node), " ", node.name)
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg
end