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

# Rule set for forward propagation
forwardAdditionMRule{T<:Number}(m_x::Array{T, 1}, m_y::Array{T, 1}) = m_x + m_y
forwardAdditionVRule{T<:Number}(V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x + V_y
forwardAdditionWRule{T<:Number}(W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x * pinv(W_x + W_y) * W_y
forwardAdditionXiRule{T<:Number}(V_x::Array{T, 2}, xi_x::Array{T, 1}, V_y::Array{T, 2}, xi_y::Array{T, 1}) = pinv(V_x + V_y) * (V_x*xi_x + V_y*xi_y)

# Calculations for a gaussian message type; Korl (2005), table 4.1
function calculateMessage!( outbound_interface_id::Int,
                            node::AdditionNode,
                            inbound_messages::Array{GaussianMessage,1})
    if outbound_interface_id == 3
        # Forward message
        msg_in_1 = inbound_messages[1]
        msg_in_2 = inbound_messages[2]
        msg_1 = deepcopy(msg_in_1)
        msg_2 = deepcopy(msg_in_2)
        msg = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1.m != nothing && msg_1.V != nothing && msg_2.m != nothing && msg_2.V != nothing
            msg.m = forwardAdditionMRule(msg_1.m, msg_2.m)
            msg.V = forwardAdditionVRule(msg_1.V, msg_2.V)
            msg.W = nothing
            msg.xi = nothing
        elseif msg_1.m != nothing && msg_1.W != nothing && msg_2.m != nothing && msg_2.W != nothing
            msg.m = forwardAdditionMRule(msg_1.m, msg_2.m)
            msg.V = nothing
            msg.W = forwardAdditionWRule(msg_1.W, msg_2.W)
            msg.xi = nothing
        elseif msg_1.xi != nothing && msg_1.V != nothing && msg_2.xi != nothing && msg_2.V != nothing
            msg.m = nothing
            msg.V = forwardAdditionVRule(msg_1.V, msg_2.V)
            msg.W = nothing
            msg.xi = forwardAdditionXiRule(msg_1.V, msg_1.xi, msg_2.V, msg_2.xi)
        #elseif msg_1.xi != nothing && msg_1.W != nothing && msg_2.xi != nothing && msg_2.W != nothing
        #    msg.m = nothing
        #    msg.V = forwardAdditionVRule(inv(msg_1.W), inv(msg_2.W))
        #    msg.W = nothing
        #    msg.xi = forwardAdditionXiRule(inv(msg_1.W), msg_1.xi, inv(msg_2.W), msg_2.xi)
        else
            error("Insufficient input to calculate outbound message on interface ", outbound_interface_id, " of ", typeof(node), " ", node.name)
        end
    elseif outbound_interface_id == 1 || outbound_interface_id == 2
        # Backward message
    end

    return node.interfaces[outbound_interface_id].message = msg
end

# ############################################
# # GeneralMessage methods
# ############################################

# Calculations for a general message type
function calculateMessage!( outbound_interface_id::Int,
                            node::AdditionNode,
                            inbound_messages::Array{GeneralMessage,1})
    if outbound_interface_id == 1 || outbound_interface_id == 2
        # Backward message
    elseif outbound_interface_id == 3
        # Forward message
        msg = GeneralMessage(inbound_messages[1].value + inbound_messages[2].value)
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg
end