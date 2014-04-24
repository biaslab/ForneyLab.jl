############################################
# FixedGainNode
############################################
# Description:
#   Multiplication with a predefined gain matrix A:
#   out = A * in1
#   Example:
#       MatrixMulitiplicationNode([1.0]; name="my_node")
#   Gain A may only be defined through the constructor!
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       GaussianMessage
#       GeneralMessage
#   2. (out):
#       GaussianMessage
#       GeneralMessage
############################################

export FixedGainNode

type FixedGainNode <: Node
    A::Array
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    out::Interface
    function FixedGainNode(A::Array; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
        self = new(ensureMatrix(deepcopy(A)), name, Array(Interface, 2))
        if !isposdef(self.A)
            warn("The specified multiplier for ", typeof(self), " ", self.name, " is not invertible. This might cause problems. Please check if this is really what you want.")
        end
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.out = self.interfaces[2]
        return self
    end
end
FixedGainNode(; args...) = FixedGainNode([1.0]; args...)

############################################
# GaussianMessage methods
############################################

# Rule set for backward propagation
backwardFixedGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}) = inv(A) * m
backwardFixedGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}, W::Array{T, 2}) = pinv(A' * W * A) * A' * W * m
backwardFixedGainVRule{T<:Number}(A::Array{T, 2}, V::Array{T, 2}) = inv(A) * V * inv(A)' # TODO: inversion only once
backwardFixedGainWRule{T<:Number}(A::Array{T, 2}, W::Array{T, 2}) = A' * W * A
backwardFixedGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}) = A' * xi

# Rule set for forward propagation
forwardFixedGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}) = A * m
forwardFixedGainVRule{T<:Number}(A::Array{T, 2}, V::Array{T, 2}) = A * V * A'
forwardFixedGainWRule{T<:Number}(A::Array{T, 2}, W::Array{T, 2}) = inv(A)' * W * inv(A) # TODO: inversion only once
forwardFixedGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}) = inv(A)' * xi
forwardFixedGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}, V::Array{T, 2}) = pinv(A * V * A') * A * V * xi # Combination of xi and V

# Calculations for a gaussian message type; Korl (2005), table 4.1
function calculateMessage!( outbound_interface_id::Int,
                            node::FixedGainNode,
                            inbound_messages::Array{GaussianMessage,1})
    if outbound_interface_id == 1
        # Backward message
        msg_in = inbound_messages[2]
        msg = deepcopy(msg_in)

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg.xi != nothing && msg.W != nothing
            msg.m = nothing
            msg.V = nothing
            msg.W = backwardFixedGainWRule(node.A, msg.W)
            msg.xi = backwardFixedGainXiRule(node.A, msg.xi)
        elseif msg.m != nothing && msg.W != nothing && isposdef(node.A) && isposdef(msg.W)
            msg.m = backwardFixedGainMRule(node.A, msg.m) # Short version of the rule
            msg.V = nothing
            msg.W = backwardFixedGainWRule(node.A, msg.W)
            msg.xi = nothing
        elseif msg.m != nothing && msg.W != nothing
            msg.m = backwardFixedGainMRule(node.A, msg.m, msg.W) # Long version of the rule
            msg.V = nothing
            msg.W = backwardFixedGainWRule(node.A, msg.W)
            msg.xi = nothing
        elseif msg.m != nothing && msg.V != nothing && isposdef(node.A) && (msg.V==zeros(size(msg.V)) || isposdef(msg.V)) # V positive definite <=> W (its inverse) is also positive definite. Border-case is an all-zero V, in which case the outbound message also has zero variance.
            msg.m = backwardFixedGainMRule(node.A, msg.m) # Short version of the rule, only valid if A is positive definite and (W is positive definite or V is 0)
            msg.V = backwardFixedGainVRule(node.A, msg.V)
            msg.W = nothing
            msg.xi = nothing
        else
            error("Insufficient input to calculate outbound message on interface ", outbound_interface_id, " of ", typeof(node), " ", node.name)
        end
    elseif outbound_interface_id == 2
        # Forward message
        msg_in = inbound_messages[1]
        msg = deepcopy(msg_in)

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg.m != nothing && msg.V != nothing
            msg.m = forwardFixedGainMRule(node.A, msg.m)
            msg.V = forwardFixedGainVRule(node.A, msg.V)
            msg.W = nothing
            msg.xi = nothing
        elseif msg.m != nothing && msg.W != nothing && isposdef(node.A)
            msg.m = forwardFixedGainMRule(node.A, msg.m)
            msg.V = nothing
            msg.W = forwardFixedGainWRule(node.A, msg.W)
            msg.xi = nothing
        elseif msg.xi != nothing && msg.W != nothing && isposdef(node.A) && isposdef(msg.W) # V should be positive definite, meaning V is invertible and its inverse W is also positive definite.
            msg.m = nothing
            msg.V = nothing
            msg.W = forwardFixedGainWRule(node.A, msg.W)
            msg.xi = forwardFixedGainXiRule(node.A, msg.xi) # Short version of the rule, only valid if A and V are positive definite
        elseif msg.xi != nothing && msg.V != nothing && msg.W != nothing && isposdef(node.A)
            msg.m = nothing
            msg.V = nothing
            msg.W = forwardFixedGainWRule(node.A, msg.W)
            msg.xi = forwardFixedGainXiRule(node.A, msg.xi, msg.V) # Long version of the rule
        else
            error("Insufficient input to calculate outbound message on interface ", outbound_interface_id, " of ", typeof(node), " ", node.name)
        end
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg
end

############################################
# GeneralMessage methods
############################################

# Calculations for a general message type
function calculateMessage!( outbound_interface_id::Int,
                            node::FixedGainNode,
                            inbound_messages::Array{GeneralMessage,1})
    if outbound_interface_id == 1
        # Backward message
        msg = GeneralMessage(pinv(node.A) * inbound_messages[2].value)
    elseif outbound_interface_id == 2
        # Forward message
        msg = GeneralMessage(node.A * inbound_messages[1].value)
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg
end