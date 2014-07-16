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
#       Message{Float64}
#       Message{Array{Float64}}
#       Message{GaussianDistribution}
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   2. (none):
#       Message{Float64}
#       Message{Array{Float64}}
#       Message{GaussianDistribution}
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   3. (none):
#       Message{Float64}
#       Message{Array{Float64}}
#       Message{GaussianDistribution}
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   ...
#   N. (none):
#       Message{Float64}
#       Message{Array{Float64}}
#       Message{GaussianDistribution}
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
############################################

export EqualityNode

type EqualityNode <: Node
    num_interfaces::Uint16
    interfaces::Array{Interface,1}
    name::ASCIIString

    function EqualityNode(num_interfaces::Integer=3; name="unnamed")
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
            return interface_id
        end
    end
    error("No free interface on $(typeof(node)) $(node.name)")
end

############################################
# GaussianDistribution methods
############################################

# Rule set
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
equalityMRule{T<:Number}(m_x::Array{T, 1}, m_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = pinv(W_x+W_y)*(W_x*m_x+W_y*m_y)
equalityVRule{T<:Number}(V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x*pinv(V_x+V_y)*V_y
equalityWRule{T<:Number}(W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x+W_y
equalityXiRule{T<:Number}(xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x+xi_y

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_value_types::Type{GaussianDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).value

    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1

    if length(node.interfaces)>3
        # Always calculate the (xi,W) parametrization of all incoming messages if it's not already present.

        # Create accumulators for W and xi of the correct size
        if !is(node.interfaces[first_incoming_id].partner.message.value.W, nothing)
            W_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.value.W))
        else
            W_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.value.V))
        end
        if !is(node.interfaces[first_incoming_id].partner.message.value.xi, nothing)
            xi_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.value.xi))
        else
            xi_sum = zeros(size(node.interfaces[first_incoming_id].partner.message.value.m))
        end

        # Sum W and xi of all incoming messages
        for incoming_interface_id = 1:length(node.interfaces)
            if incoming_interface_id==outbound_interface_id
                continue
            end
            # Calculate (xi,W) parametrization if it's not available yet
            ensureXiWParametrization!(node.interfaces[incoming_interface_id].partner.message.value)
            W_sum += node.interfaces[incoming_interface_id].partner.message.value.W
            xi_sum += node.interfaces[incoming_interface_id].partner.message.value.xi
        end
        dist_out.xi = xi_sum
        dist_out.W = W_sum
        dist_out.V = nothing
        dist_out.m = nothing
    else # 3 interfaces
        dist_1 = node.interfaces[first_incoming_id].partner.message.value
        dist_2 = node.interfaces[(outbound_interface_id==3) ? 2 : 3].partner.message.value
        if dist_1.m != nothing && dist_1.W != nothing && dist_2.m != nothing && dist_2.W != nothing
            dist_out.m  = equalityMRule(dist_1.m, dist_2.m, dist_1.W, dist_2.W)
            dist_out.V  = nothing
            dist_out.W  = equalityWRule(dist_1.W, dist_2.W)
            dist_out.xi = nothing
        elseif dist_1.xi != nothing && dist_1.V != nothing && dist_2.xi != nothing && dist_2.V != nothing
            dist_out.m  = nothing
            dist_out.V  = equalityVRule(dist_1.V, dist_2.V)
            dist_out.W  = nothing
            dist_out.xi = equalityXiRule(dist_1.xi, dist_2.xi)
        else
            # Use (xi,W)
            ensureXiWParametrization!(dist_1)
            ensureXiWParametrization!(dist_2)
            dist_out.m  = nothing
            dist_out.V  = nothing
            dist_out.W  = equalityWRule(dist_1.W, dist_2.W)
            dist_out.xi = equalityXiRule(dist_1.xi, dist_2.xi)
        end
    end

    return node.interfaces[outbound_interface_id].message
end

############################################
# Float64 and Array{Float64} methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_value_types::Union(Type{Float64}, Type{Array{Float64}}, Type{Array{Float64, 2}}))
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    # Outbound message is equal to the inbound messages if not all inbound messages are equal.
    # Otherwise, the outbound message is Message{Float64}(0.0)
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    incoming_msg = node.interfaces[first_incoming_id].partner.message
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id || interface_id==first_incoming_id
            continue
        end
        if node.interfaces[interface_id].partner.message.value != node.interfaces[first_incoming_id].partner.message.value
            if typeof(node.interfaces[first_incoming_id].partner.message.value)<:Array
                ans = zeros(size(incoming_msg.value))
            else
                ans = 0.0
            end
            msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], inbound_messages_value_types, size(ans))
            msg_out.value = ans
            return node.interfaces[outbound_interface_id].message
        end
    end
    ans = deepcopy(incoming_msg.value)
    msg_out = getOrCreateMessage(node.interfaces[outbound_interface_id], inbound_messages_value_types, size(ans))
    msg_out.value = ans

    return node.interfaces[outbound_interface_id].message
end

############################################
# InverseGammaDistribution methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_value_types::Type{InverseGammaDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], InverseGammaDistribution).value

    # Definition from Korl table 5.2
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    if length(node.interfaces)!=3 
        error("Equality update rule for inverse gamma distribution only defined for three interfaces")
    end
    dist_out.a = 1.0
    dist_out.b = 0.0
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id
            continue
        end
        dist_out.a += node.interfaces[interface_id].partner.message.value.a
        dist_out.b += node.interfaces[interface_id].partner.message.value.b
    end

    return node.interfaces[outbound_interface_id].message
end

############################################
# GammaDistribution methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages_value_types::Type{GammaDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GammaDistribution).value

    # Definition from derivation in notebook gamma_message_eq_node_derivation
    first_incoming_id = (outbound_interface_id==1) ? 2 : 1
    if length(node.interfaces)!=3 
        error("Equality update rule for gamma distribution only defined for three interfaces")
    end
    dist_out.a = -1.0
    dist_out.b = 0.0
    for interface_id = 1:length(node.interfaces)
        if interface_id==outbound_interface_id
            continue
        end
        dist_out.a += node.interfaces[interface_id].partner.message.value.a
        dist_out.b += node.interfaces[interface_id].partner.message.value.b
    end

    return node.interfaces[outbound_interface_id].message
end
