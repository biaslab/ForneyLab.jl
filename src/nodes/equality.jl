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
#   2. (none):
#       GaussianMessage
#       GeneralMessage
#   3. (none):
#       GaussianMessage
#       GeneralMessage
#   ...
#   N. (none):
#       GaussianMessage
#       GeneralMessage
############################################

export EqualityNode

type EqualityNode <: Node
    num_interfaces::Uint16
    interfaces::Array{Interface,1}
    name::ASCIIString

    function EqualityNode(num_interfaces::Integer; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        @assert(num_interfaces>2, "An EqualityNode should have at least 3 interfaces")
        self = new(num_interfaces, Array(Interface, num_interfaces), name)
        # Create interfaces
        for interface_id=1:num_interfaces
            self.interfaces[interface_id] = Interface(self)
        end
        return self
    end
end
EqualityNode(; args...) = EqualityNode(3; args...)

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
                            inbound_messages::Array{GaussianMessage, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end

    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1

    if length(node.interfaces)>3
        # Always calculate the (xi,W) parametrization of all incoming messages if it's not already present.

        # Create accumulators for W and xi of the correct size
        if !is(inbound_messages[first_incoming_id].W, nothing)
            W_sum = zeros(size(inbound_messages[first_incoming_id].W))
        else
            W_sum = zeros(size(inbound_messages[first_incoming_id].V))
        end
        if !is(inbound_messages[first_incoming_id].xi, nothing)
            xi_sum = zeros(size(inbound_messages[first_incoming_id].xi))
        else
            xi_sum = zeros(size(inbound_messages[first_incoming_id].m))
        end

        # Sum W and xi of all incoming messages
        for incoming_interface_id = 1:length(node.interfaces)
            if incoming_interface_id==outbound_interface_id
                continue
            end
            # Calculate (xi,W) parametrization if it's not available yet
            ensureXiWParametrization!(inbound_messages[incoming_interface_id])
            W_sum += inbound_messages[incoming_interface_id].W
            xi_sum += inbound_messages[incoming_interface_id].xi
        end
        msg_out = GaussianMessage(xi=xi_sum, W=W_sum)
    else # 3 interfaces
        msg_1 = inbound_messages[first_incoming_id]
        msg_2 = inbound_messages[(outbound_interface_id==3) ? 2 : 3]
        msg_out = GaussianMessage()
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

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg_out
end

############################################
# GeneralMessage methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages::Array{GeneralMessage, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end

    # Outbound message is equal to the inbound messages if not all inbound messages are equal.
    # Otherwise, the outbound message is GeneralMessage(0.0)
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id || interface_id==first_incoming_id
            continue
        end
        if inbound_messages[interface_id].value!=inbound_messages[first_incoming_id].value
            if typeof(inbound_messages[first_incoming_id].value)<:Array
                return node.interfaces[outbound_interface_id].message = GeneralMessage(zeros(size(inbound_messages[first_incoming_id].value)))
            else
                return node.interfaces[outbound_interface_id].message = GeneralMessage(0.0)
            end
        end
    end

    return node.interfaces[outbound_interface_id].message = deepcopy(inbound_messages[first_incoming_id])
end