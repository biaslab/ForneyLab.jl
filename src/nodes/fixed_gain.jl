############################################
# FixedGainNode
############################################
# Description:
#   Multiplication: out = A*in
#
#     in      out
#   ----->[A]----->
#
#   f(in,out) = δ(out - A*in)
#
# Interfaces:
#   1 i[:in], 2 i[:out]
#
# Construction:
#   FixedGainNode([1.0], name="my_node")
#
############################################

export FixedGainNode

type FixedGainNode <: Node
    A::Array
    name::ASCIIString
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible

    function FixedGainNode(A::Union(Array{Float64},Float64)=1.0; name=unnamedStr())
        # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
        A = (typeof(A)==Float64) ? fill!(Array(Float64,1,1),A) : ensureMatrix(deepcopy(A))
        self = new(A, name, Array(Interface, 2), Dict{Symbol,Interface}())

        for (iface_id, iface_name) in enumerate([:in, :out])
            self.i[iface_name] = self.interfaces[iface_id] = Interface(self)
        end

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

# Backward Gaussian to IN
function fixedGainGaussianBackwardRule!(dist_result::GaussianDistribution, dist_2::GaussianDistribution, A::Any, A_inv::Any=nothing)
    # Calculations for a gaussian message type; Korl (2005), table 4.1
    # dist_result = inv(A) * dist_2

    if isValid(dist_2.xi) && isValid(dist_2.W)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardFixedGainWRule(A, dist_2.W)
        dist_result.xi = backwardFixedGainXiRule(A, dist_2.xi)
    elseif isValid(dist_2.m) && isValid(dist_2.W) && isRoundedPosDef(A) && isRoundedPosDef(dist_2.W)
        dist_result.m = backwardFixedGainMRule(A_inv, dist_2.m) # Short version of the rule
        invalidate!(dist_result.V)
        dist_result.W = backwardFixedGainWRule(A, dist_2.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_2.m) && isValid(dist_2.W)
        dist_result.m = backwardFixedGainMRule(A, dist_2.m, dist_2.W) # Long version of the rule
        invalidate!(dist_result.V)
        dist_result.W = backwardFixedGainWRule(A, dist_2.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_2.m) && isValid(dist_2.V) && (A_inv != nothing) && (dist_2.V==zeros(size(dist_2.V)) || isRoundedPosDef(dist_2.V)) # V positive definite <=> W (its inverse) is also positive definite. Border-case is an all-zero V, in which case the outbound message also has zero variance.
        dist_result.m = backwardFixedGainMRule(A_inv, dist_2.m) # Short version of the rule, only valid if A is positive definite and (W is positive definite or V is 0)
        dist_result.V = backwardFixedGainVRule(A_inv, dist_2.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    else
        # Fallback: convert inbound message to (xi,W) parametrization and then use efficient rules
        ensureXiWParametrization!(dist_2)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardFixedGainWRule(A, dist_2.W)
        dist_result.xi = backwardFixedGainXiRule(A, dist_2.xi)
    end   

    return dist_result 
end

function sumProduct!(node::FixedGainNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            msg_out::Message{GaussianDistribution})

    dist_1 = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if outbound_interface_id == 1
        fixedGainGaussianBackwardRule!(dist_1, msg_out.payload, node.A, isdefined(node, :A_inv) ? node.A_inv : nothing)
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end
    
    return (:fixed_gain_gaussian_backward,
            node.interfaces[outbound_interface_id].message)
end

# Forward Gaussian to OUT
function fixedGainGaussianForwardRule!(dist_result::GaussianDistribution, dist_1::GaussianDistribution, A::Any, A_inv::Any=nothing)
    # Calculations for a gaussian message type; Korl (2005), table 4.1
    # dist_result = A * dist_1

    if isValid(dist_1.m) && isValid(dist_1.V)
        dist_result.m = forwardFixedGainMRule(A, dist_1.m)
        dist_result.V = forwardFixedGainVRule(A, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.m) && isValid(dist_1.W) && (A_inv != nothing)
        dist_result.m = forwardFixedGainMRule(A, dist_1.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardFixedGainWRule(A_inv, dist_1.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.xi) && isValid(dist_1.W) && isRoundedPosDef(A) && isRoundedPosDef(dist_1.W) && (A_inv != nothing) # V should be positive definite, meaning V is invertible and its inverse W is also positive definite.
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardFixedGainWRule(A_inv, dist_1.W)
        dist_result.xi = forwardFixedGainXiRule(A_inv, dist_1.xi) # Short version of the rule, only valid if A and V are positive definite
    elseif isValid(dist_1.xi) && isValid(dist_1.V) && isValid(dist_1.W) && (A_inv != nothing)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardFixedGainWRule(A_inv, dist_1.W)
        dist_result.xi = forwardFixedGainXiRule(A_inv, dist_1.xi, dist_1.V) # Long version of the rule
    else
        # Fallback: convert inbound message to (m,V) parametrization and then use efficient rules
        ensureMVParametrization!(dist_1)
        dist_result.m = forwardFixedGainMRule(A, dist_1.m)
        dist_result.V = forwardFixedGainVRule(A, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    end

    return dist_result 
end

function sumProduct!(node::FixedGainNode,
                            outbound_interface_id::Int,
                            msg_in::Message{GaussianDistribution},
                            ::Nothing)
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if outbound_interface_id == 2
        # Forward message
        dist_1 = msg_in.payload

        fixedGainGaussianForwardRule!(dist_2, msg_in.payload, node.A, isdefined(node, :A_inv) ? node.A_inv : nothing)
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end
    
    return (:fixed_gain_gaussian_forward,
            node.interfaces[outbound_interface_id].message)
end

# Backward DeltaDistribution to IN
function sumProduct!{T<:Any}(
                            node::FixedGainNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            msg_out::Message{DeltaDistribution{T}})
    if outbound_interface_id == 1
        # Backward message
        ans = node.A_inv * msg_out.payload.m
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end

    msg_ans = ensureMessage!(node.interfaces[outbound_interface_id], DeltaDistribution{typeof(ans)})
    msg_ans.payload.m = ans

    return (:fixed_gain_delta_backward,
            node.interfaces[outbound_interface_id].message)
end

# Forward DeltaDistribution to OUT
function sumProduct!{T<:Any}(
                            node::FixedGainNode,
                            outbound_interface_id::Int,
                            msg_in::Message{DeltaDistribution{T}},
                            ::Nothing)
    if outbound_interface_id == 2
        # Forward message
        ans = node.A * msg_in.payload.m
    else
        error("Invalid interface id $(outbound_interface_id) for calculating message on $(typeof(node)) $(node.name)")
    end

    msg_ans = ensureMessage!(node.interfaces[outbound_interface_id], DeltaDistribution{typeof(ans)})
    msg_ans.payload.m = ans

    return (:fixed_gain_delta_forward,
            node.interfaces[outbound_interface_id].message)
end