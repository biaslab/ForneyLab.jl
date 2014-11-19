############################################
# EqualityNode
############################################
# Description:
#   Equality constraint node with 3 (symmetrical) interfaces.
#
#   Example:
#       EqualityNode(name="my_equ")
#
# Interface ids, (names) and supported message types:
# 1,2,3. (none):
#       Message{DeltaDistribution}
#       Message{GaussianDistribution}
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#       Message{NormalGammaDistribution}
#       Message{StudentsTDistribution}
############################################

export EqualityNode

type EqualityNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}

    function EqualityNode(; name=unnamedStr())
        self = new(name, Array(Interface, 3))
        # Create interfaces
        for interface_id = 1:3
            self.interfaces[interface_id] = Interface(self)
        end
        return self
    end
end

isDeterministic(::EqualityNode) = true

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

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_1::Message{GaussianDistribution},
                            msg_2::Message{GaussianDistribution},
                            ::Nothing)
    # Calculate an outbound message based on the inbound messages and the node function.
    dist_result = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
    
    dist_1 = msg_1.payload
    dist_2 = msg_2.payload
    if dist_1.m != nothing && dist_1.W != nothing && dist_2.m != nothing && dist_2.W != nothing
        dist_result.m  = equalityMRule(dist_1.m, dist_2.m, dist_1.W, dist_2.W)
        dist_result.V  = nothing
        dist_result.W  = equalityWRule(dist_1.W, dist_2.W)
        dist_result.xi = nothing
    elseif dist_1.xi != nothing && dist_1.V != nothing && dist_2.xi != nothing && dist_2.V != nothing
        dist_result.m  = nothing
        dist_result.V  = equalityVRule(dist_1.V, dist_2.V)
        dist_result.W  = nothing
        dist_result.xi = equalityXiRule(dist_1.xi, dist_2.xi)
    else
        # Use (xi,W)
        ensureXiWParametrization!(dist_1)
        ensureXiWParametrization!(dist_2)
        dist_result.m  = nothing
        dist_result.V  = nothing
        dist_result.W  = equalityWRule(dist_1.W, dist_2.W)
        dist_result.xi = equalityXiRule(dist_1.xi, dist_2.xi)
    end

    return node.interfaces[outbound_interface_id].message
end
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_1::Message{GaussianDistribution}, msg_2::Message{GaussianDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_1::Message{GaussianDistribution}, ::Nothing, msg_2::Message{GaussianDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_st::Message{StudentsTDistribution},
                            msg_n::Message{GaussianDistribution},
                            ::Nothing)
    # Combination of Gaussian and student's t-distribution
    # Definitions available in derivations notebook
    # Same as Gaussian equality rule

    dist_result = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    dist_1 = GaussianDistribution(m=msg_st.payload.m, W=msg_st.payload.W) # Approximate student's with Gaussian
    dist_2 = msg_n.payload

    if dist_1.m != nothing && dist_1.W != nothing && dist_2.m != nothing && dist_2.W != nothing
        dist_result.m  = equalityMRule(dist_1.m, dist_2.m, dist_1.W, dist_2.W)
        dist_result.V  = nothing
        dist_result.W  = equalityWRule(dist_1.W, dist_2.W)
        dist_result.xi = nothing
    elseif dist_1.xi != nothing && dist_1.V != nothing && dist_2.xi != nothing && dist_2.V != nothing
        dist_result.m  = nothing
        dist_result.V  = equalityVRule(dist_1.V, dist_2.V)
        dist_result.W  = nothing
        dist_result.xi = equalityXiRule(dist_1.xi, dist_2.xi)
    else
        # Use (xi,W)
        ensureXiWParametrization!(dist_1)
        ensureXiWParametrization!(dist_2)
        dist_result.m  = nothing
        dist_result.V  = nothing
        dist_result.W  = equalityWRule(dist_1.W, dist_2.W)
        dist_result.xi = equalityXiRule(dist_1.xi, dist_2.xi)
    end

    return node.interfaces[outbound_interface_id].message
end
# Call signature for messages other ways around
# Outbound interface 3
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_n::Message{GaussianDistribution}, msg_st::Message{StudentsTDistribution}, msg_dummy::Nothing) = updateNodeMessage!(node, outbound_interface_id, msg_st, msg_n, msg_dummy)
# Outbound interface 2
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_n::Message{GaussianDistribution}, msg_dummy::Nothing, msg_st::Message{StudentsTDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_st, msg_n, msg_dummy)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_st::Message{StudentsTDistribution}, msg_dummy::Nothing, msg_n::Message{GaussianDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_st, msg_n, msg_dummy)
# Outbound interface 1
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_dummy::Nothing, msg_n::Message{GaussianDistribution}, msg_st::Message{StudentsTDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_st, msg_n, msg_dummy)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_dummy::Nothing, msg_st::Message{StudentsTDistribution}, msg_n::Message{GaussianDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_st, msg_n, msg_dummy)

############################################
# DeltaDistribution methods
############################################

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_1::Message{DeltaDistribution},
                            msg_2::Message{DeltaDistribution},
                            ::Nothing)
    # Calculate an outbound message based on the inbound messages and the node function.

    # Outbound message is equal to the inbound messages if both inbound messages are equal.
    # Otherwise, the outbound message is 0.0

    msg_result = getOrCreateMessage(node.interfaces[outbound_interface_id], DeltaDistribution)
    if msg_1.payload == msg_2.payload
        msg_result.payload = copy(msg_1.payload)
    else
        msg_result.payload = zero(msg_1.payload)
    end

    return node.interfaces[outbound_interface_id].message
end
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_1::Message{DeltaDistribution}, ::Nothing, msg_2::Message{DeltaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_1::Message{DeltaDistribution}, msg_2::Message{DeltaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)


############################################
# InverseGammaDistribution methods
############################################

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_1::Message{InverseGammaDistribution},
                            msg_2::Message{InverseGammaDistribution},
                            ::Nothing)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload

    # Definition from Korl table 5.2
    dist_out.a = 1.0 + msg_1.payload.a + msg_2.payload.a
    dist_out.b = msg_1.payload.b + msg_2.payload.b

    return node.interfaces[outbound_interface_id].message
end
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_1::Message{InverseGammaDistribution}, ::Nothing, msg_2::Message{InverseGammaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_1::Message{InverseGammaDistribution}, msg_2::Message{InverseGammaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)

############################################
# GammaDistribution methods
############################################

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_1::Message{GammaDistribution},
                            msg_2::Message{GammaDistribution},
                            ::Nothing)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GammaDistribution).payload

    # Derivation available in notebook
    dist_out.a = -1.0 + msg_1.payload.a + msg_2.payload.a
    dist_out.b = msg_1.payload.b + msg_2.payload.b

    return node.interfaces[outbound_interface_id].message
end
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_1::Message{GammaDistribution}, ::Nothing, msg_2::Message{GammaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_1::Message{GammaDistribution}, msg_2::Message{GammaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_1, msg_2, nothing)

############################################
# Gaussian-DeltaDistribution combination
############################################

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_delta::Message{DeltaDistribution},
                            msg_n::Message{GaussianDistribution},
                            ::Nothing)
    # Combination of Gaussian and float
    return node.interfaces[outbound_interface_id].message = deepcopy(msg_delta)
end
# Call signature for messages other ways around
# Outbound interface 3
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_n::Message{GaussianDistribution}, msg_delta::Message{DeltaDistribution}, ::Nothing) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_n, nothing)
# Outbound interface 2
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_n::Message{GaussianDistribution}, ::Nothing, msg_delta::Message{DeltaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_n, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_delta::Message{DeltaDistribution}, ::Nothing, msg_n::Message{GaussianDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_n, nothing)
# Outbound interface 1
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_n::Message{GaussianDistribution}, msg_delta::Message{DeltaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_n, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_delta::Message{DeltaDistribution}, msg_n::Message{GaussianDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_n, nothing)

############################################
# Gamma-DeltaDistribution combination
############################################

function updateNodeMessage!(node::EqualityNode,
                            outbound_interface_id::Int,
                            msg_delta::Message{DeltaDistribution},
                            msg_gam::Message{GammaDistribution},
                            ::Nothing)
    # Combination of gamma and float
    (typeof(msg_delta.payload.m) <: Number) || error("Node $(node.name) cannot perform update for gamma-delta combination for non-scalars")
    (msg_delta.payload.m >= 0) || error("Node $(node.name) cannot perform update for gamma-delta combination for negative numbers")
    return node.interfaces[outbound_interface_id].message = deepcopy(msg_delta)
end
# Call signature for messages other ways around
# Outbound interface 3
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_gam::Message{GammaDistribution}, msg_delta::Message{DeltaDistribution}, ::Nothing) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_gam, nothing)
# Outbound interface 2
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_gam::Message{GammaDistribution}, ::Nothing, msg_delta::Message{DeltaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_gam, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, msg_delta::Message{DeltaDistribution}, ::Nothing, msg_gam::Message{GammaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_gam, nothing)
# Outbound interface 1
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_gam::Message{GammaDistribution}, msg_delta::Message{DeltaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_gam, nothing)
updateNodeMessage!(node::EqualityNode, outbound_interface_id::Int, ::Nothing, msg_delta::Message{DeltaDistribution}, msg_gam::Message{GammaDistribution}) = updateNodeMessage!(node, outbound_interface_id, msg_delta, msg_gam, nothing)