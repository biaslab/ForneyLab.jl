############################################
# FixedGainNode
############################################
# Description:
#   Multiplication with a predefined gain matrix A:
#
#    in1      out
#   ----->[A]----->
#
#   out = A * in1
#
#   Example:
#       FixedGainNode([1.0]; name="my_node")
#   Gain A may only be defined through the constructor!
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       Message{DeltaDistribution}
#       Message{GaussianDistribution}
#   2. (out):
#       Message{DeltaDistribution}
#       Message{GaussianDistribution}
############################################

export FixedGainNode

type FixedGainNode <: Node
    A::Array
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Helper fields filled by constructor
    in1::Interface
    out::Interface
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible
    function FixedGainNode(A::Union(Array{Float64},Float64)=1.0; name=unnamedStr())
        # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
        A = (typeof(A)==Float64) ? fill!(Array(Float64,1,1),A) : ensureMatrix(deepcopy(A))
        self = new(A, name, Array(Interface, 2))

        # Set up the interfaces
        self.in1 = self.interfaces[1] = Interface(self)
        self.out = self.interfaces[2] = Interface(self)

        # Try to precompute inv(A)
        try
            self.A_inv = inv(self.A)
        catch
            warn("The specified multiplier for $(typeof(self)) $(self.name) is not invertible. This might cause problems. Double check that this is what you want.")
        end
        return self
    end
end

isDeterministic(::FixedGainNode) = true

############################################
# GaussianDistribution methods
############################################

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardFixedGainMRule{T<:Number}(A_inv::Array{T, 2}, m::Array{T, 1}) = A_inv * m
backwardFixedGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}, W::Array{T, 2}) = pinv(A' * W * A) * A' * W * m
backwardFixedGainVRule{T<:Number}(A_inv::Array{T, 2}, V::Array{T, 2}) = A_inv * V * A_inv'
backwardFixedGainWRule{T<:Number}(A::Array{T, 2}, W::Array{T, 2}) = A' * W * A
backwardFixedGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}) = A' * xi

# Rule set for forward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardFixedGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}) = A * m
forwardFixedGainVRule{T<:Number}(A::Array{T, 2}, V::Array{T, 2}) = A * V * A'
forwardFixedGainWRule{T<:Number}(A_inv::Array{T, 2}, W::Array{T, 2}) = A_inv' * W * A_inv
forwardFixedGainXiRule{T<:Number}(A_inv::Array{T, 2}, xi::Array{T, 1}) = A_inv' * xi
forwardFixedGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}, V::Array{T, 2}) = pinv(A * V * A') * A * V * xi # Combination of xi and V

# Backward Gaussian to IN1
function updateNodeMessage!(node::FixedGainNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            msg_out::Message{GaussianDistribution})

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    # Calculations for a gaussian message type; Korl (2005), table 4.1
    if outbound_interface_id == 1
        dist_2 = msg_out.payload
        # Select parameterization
        # Order is from least to most computationally intensive
        if dist_2.xi != nothing && dist_2.W != nothing
            dist_out.m = nothing
            dist_out.V = nothing
            dist_out.W = backwardFixedGainWRule(node.A, dist_2.W)
            dist_out.xi = backwardFixedGainXiRule(node.A, dist_2.xi)
        elseif dist_2.m != nothing && dist_2.W != nothing && isRoundedPosDef(node.A) && isRoundedPosDef(dist_2.W)
            dist_out.m = backwardFixedGainMRule(node.A_inv, dist_2.m) # Short version of the rule
            dist_out.V = nothing
            dist_out.W = backwardFixedGainWRule(node.A, dist_2.W)
            dist_out.xi = nothing
        elseif dist_2.m != nothing && dist_2.W != nothing
            dist_out.m = backwardFixedGainMRule(node.A, dist_2.m, dist_2.W) # Long version of the rule
            dist_out.V = nothing
            dist_out.W = backwardFixedGainWRule(node.A, dist_2.W)
            dist_out.xi = nothing
        elseif dist_2.m != nothing && dist_2.V != nothing && isdefined(node, :A_inv) && (dist_2.V==zeros(size(dist_2.V)) || isRoundedPosDef(dist_2.V)) # V positive definite <=> W (its inverse) is also positive definite. Border-case is an all-zero V, in which case the outbound message also has zero variance.
            dist_out.m = backwardFixedGainMRule(node.A_inv, dist_2.m) # Short version of the rule, only valid if A is positive definite and (W is positive definite or V is 0)
            dist_out.V = backwardFixedGainVRule(node.A_inv, dist_2.V)
            dist_out.W = nothing
            dist_out.xi = nothing
        else
            # Fallback: convert inbound message to (xi,W) parametrization and then use efficient rules
            ensureXiWParametrization!(dist_2)
            dist_out.m = nothing
            dist_out.V = nothing
            dist_out.W = backwardFixedGainWRule(node.A, dist_2.W)
            dist_out.xi = backwardFixedGainXiRule(node.A, dist_2.xi)
        end
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end
    
    return node.interfaces[outbound_interface_id].message
end

# Forward Gaussian to OUT
function updateNodeMessage!(node::FixedGainNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{GaussianDistribution},
                            ::Nothing)

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    # Calculations for a gaussian message type; Korl (2005), table 4.1
    if outbound_interface_id == 2
        # Forward message
        dist_1 = msg_in1.payload
        # Select parameterization
        # Order is from least to most computationally intensive
        if dist_1.m != nothing && dist_1.V != nothing
            dist_out.m = forwardFixedGainMRule(node.A, dist_1.m)
            dist_out.V = forwardFixedGainVRule(node.A, dist_1.V)
            dist_out.W = nothing
            dist_out.xi = nothing
        elseif dist_1.m != nothing && dist_1.W != nothing && isdefined(node, :A_inv)
            dist_out.m = forwardFixedGainMRule(node.A, dist_1.m)
            dist_out.V = nothing
            dist_out.W = forwardFixedGainWRule(node.A_inv, dist_1.W)
            dist_out.xi = nothing
        elseif dist_1.xi != nothing && dist_1.W != nothing && isRoundedPosDef(node.A) && isRoundedPosDef(dist_1.W) # V should be positive definite, meaning V is invertible and its inverse W is also positive definite.
            dist_out.m = nothing
            dist_out.V = nothing
            dist_out.W = forwardFixedGainWRule(node.A_inv, dist_1.W)
            dist_out.xi = forwardFixedGainXiRule(node.A_inv, dist_1.xi) # Short version of the rule, only valid if A and V are positive definite
        elseif dist_1.xi != nothing && dist_1.V != nothing && dist_1.W != nothing && isdefined(node, :A_inv)
            dist_out.m = nothing
            dist_out.V = nothing
            dist_out.W = forwardFixedGainWRule(node.A_inv, dist_1.W)
            dist_out.xi = forwardFixedGainXiRule(node.A_inv, dist_1.xi, dist_1.V) # Long version of the rule
        else
            # Fallback: convert inbound message to (m,V) parametrization and then use efficient rules
            ensureMVParametrization!(dist_1)
            dist_out.m = forwardFixedGainMRule(node.A, dist_1.m)
            dist_out.V = forwardFixedGainVRule(node.A, dist_1.V)
            dist_out.W = nothing
            dist_out.xi = nothing
        end
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end
    
    return node.interfaces[outbound_interface_id].message
end

# Backward DeltaDistribution to IN1
function updateNodeMessage!(node::FixedGainNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            msg_out::Message{DeltaDistribution})
    if outbound_interface_id == 1
        # Backward message
        ans = node.A_inv * msg_out.payload.m
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end

    msg_ans = getOrCreateMessage(node.interfaces[outbound_interface_id], DeltaDistribution)
    msg_ans.payload.m = ans

    return node.interfaces[outbound_interface_id].message

end

# Forward DeltaDistribution to OUT
function updateNodeMessage!(node::FixedGainNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{DeltaDistribution},
                            ::Nothing)
    if outbound_interface_id == 2
        # Forward message
        ans = node.A * msg_in1.payload.m
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end

    msg_ans = getOrCreateMessage(node.interfaces[outbound_interface_id], DeltaDistribution)
    msg_ans.payload.m = ans

    return node.interfaces[outbound_interface_id].message
end