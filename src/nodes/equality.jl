############################################
# EqualityNode
############################################
# Description:
#   Equality constraint node with variable
#   number of (symmetrical) interfaces (ports).
#   Interfaces don't have names since the number
#   of interfaces is variable.
#
#   Example:
#       EqualityNode(name="3_port_equ") # 3 interfaces is default
#       EqualityNode(5; name="5_port_equ")
#
# Interface ids, (names) and supported message types:
#   1. (none):
#       GaussianMessage
#       GeneralMessage
#       GammaMessage
#   2. (none):
#       GaussianMessage
#       GeneralMessage
#       GammaMessage
#   3. (none):
#       GaussianMessage
#       GeneralMessage
#       GammaMessage
#   ...
#   N. (none):
#       GaussianMessage
#       GeneralMessage
#       GammaMessage
############################################

export EqualityNode

type EqualityNode <: Node
    num_interfaces::Uint16
    interfaces::Array{Interface,1}
    name::ASCIIString

    function EqualityNode(num_interfaces::Integer=3; name="unnamed", args...)
        @assert(num_interfaces>2, "An EqualityNode should have at least 3 interfaces")
        self = new(num_interfaces, Array(Interface, num_interfaces), name)
        # Create interfaces
        for interface_id=1:num_interfaces
            self.interfaces[interface_id] = Interface(self)
        end
        return self
    end
end

# Overload firstFreeInterface since EqualityNode is symmetrical in its interfaces
function firstFreeInterface(node::EqualityNode)
    # Return id of first free interface of a symmetrical node
    for interface_id = 1:length(node.interfaces)
        if node.interfaces[interface_id].partner == nothing
            return interrace_id
        end
    end
    error("No free interface on $(typeof(node)) $(node.name)")
end

############################################
# GaussianMessage methods
############################################

# Rule set
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
equalityMRule{T<:Number}(m_x::Array{T, 1}, m_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = pinv(W_x+W_y)*(W_x*m_x+W_y*m_y)
equalityVRule{T<:Number}(V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x*pinv(V_x+V_y)*V_y
equalityWRule{T<:Number}(W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x+W_y
equalityXiRule{T<:Number}(xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x+xi_y

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_types::Type{GaussianMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    msg_out = getOrAssign(node.interfaces[outbound_interface_id], GaussianMessage)

    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1

    if length(node.interfaces)>3
        # Always calculate the (xi,W) parametrization of all incoming messages if it's not already present.

        # Create accumulators for W and xi of the correct size
        if !is(node.interfaces[first_incoming_id].partner.message.W, nothing)
            W_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.W))
        else
            W_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.V))
        end
        if !is(node.interfaces[first_incoming_id].partner.message.xi, nothing)
            xi_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.xi))
        else
            xi_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.m))
        end

        # Sum W and xi of all incoming messages
        for incoming_interface_id = 1:length(node.interfaces)
            if incoming_interface_id==outbound_interface_id
                continue
            end
            # Calculate (xi,W) parametrization if it's not available yet
            ensureXiWParametrization!(node.interfaces[incoming_interface_id].partner.message)
            W_sum += node.interfaces[incoming_interface_id].partner.message.W
            xi_sum += node.interfaces[incoming_interface_id].partner.message.xi
        end
        msg_out.xi = xi_sum
        msg_out.W = W_sum
        msg_out.V = nothing
        msg_out.m = nothing
    else # 3 interfaces
        msg_1 = node.interfaces[first_incoming_id].partner.message
        msg_2 = node.interfaces[(outbound_interface_id==3) ? 2 : 3].partner.message
        if msg_1.m != nothing && msg_1.W != nothing && msg_2.m != nothing && msg_2.W != nothing
            msg_out.m  = equalityMRule(msg_1.m, msg_2.m, msg_1.W, msg_2.W)
            msg_out.V  = nothing
            msg_out.W  = equalityWRule(msg_1.W, msg_2.W)
            msg_out.xi = nothing
        elseif msg_1.xi != nothing && msg_1.V != nothing && msg_2.xi != nothing && msg_2.V != nothing
            msg_out.m  = nothing
            msg_out.V  = equalityVRule(msg_1.V, msg_2.V)
            msg_out.W  = nothing
            msg_out.xi = equalityXiRule(msg_1.xi, msg_2.xi)
        else
            # Use (xi,W)
            ensureXiWParametrization!(msg_1)
            ensureXiWParametrization!(msg_2)
            msg_out.m  = nothing
            msg_out.V  = nothing
            msg_out.W  = equalityWRule(msg_1.W, msg_2.W)
            msg_out.xi = equalityXiRule(msg_1.xi, msg_2.xi)
        end
    end

    return msg_out
end

############################################
# GeneralMessage methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_types::Type{GeneralMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    msg_out = getOrAssign(node.interfaces[outbound_interface_id], GeneralMessage)

    # Outbound message is equal to the inbound messages if not all inbound messages are equal.
    # Otherwise, the outbound message is GeneralMessage(0.0)
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id || interface_id==first_incoming_id
            continue
        end
        if node.interfaces[interface_id].partner.message.value!=node.interfaces[first_incoming_id].partner.message.value
            if typeof(node.interfaces[first_incoming_id].partner.message.value)<:Array
                msg_out.value = zeros(size(node.interfaces[first_incoming_id].partner.message.value))
            else
                msg_out.value = 0.0
            end
            return msg_out
        end
    end
    msg_out.value = deepcopy(node.interfaces[first_incoming_id].partner.message.value)

    return msg_out
end

############################################
# InverseGammaMessage methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_types::Type{InverseGammaMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    msg_out = getOrAssign(node.interfaces[outbound_interface_id], InverseGammaMessage)

    # Definition from Korl table 5.2
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    if length(node.interfaces)!=3 
        error("Equality update rule for inverse gamma distribution only defined for three interfaces")
    end
    a = 1.0
    b = 0
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id
            continue
        end
        msg = node.interfaces[interface_id].partner.message
        a += msg.a
        b += msg.b
    end
    msg_out.a = a
    msg_out.b = b

    return msg_out
end

############################################
# GammaMessage methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_types::Type{GammaMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    msg_out = getOrAssign(node.interfaces[outbound_interface_id], GammaMessage)

    # Definition from derivation in notebook gamma_message_eq_node_derivation
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    if length(node.interfaces)!=3 
        error("Equality update rule for gamma distribution only defined for three interfaces")
    end
    a = -1.0
    b = 0
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id
            continue
        end
        msg = node.interfaces[interface_id].partner.message
        a += msg.a
        b += msg.b
    end
    msg_out.a = a
    msg_out.b = b

    return msg_out
end
