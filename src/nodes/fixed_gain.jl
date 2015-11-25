############################################
# GainNode
############################################
# Description:
#   Multiplication:
#         If the message at the edge "gain" exists, then:
#              out = gain * in
#         Otherwise:
#              out = A*in
#
#         gain
#          |
#     in   V   out
#   ----->[A]----->
#
#
#                       | δ(out - A*in), if gain = nothing
#     f(in,out,gain) =  |
#                       | δ(out - gain*in), otherwise
# Interfaces:
#   1 i[:in], 2 i[:out], 3 i[:gain]
#
# Construction:
#    GainNode([1.0], id=:my_node)
#     or
#    GainNode(id=:my_node)
#
############################################

export GainNode

type GainNode <: Node
    A::Matrix{Float64}
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    A_inv::Matrix{Float64} # holds pre-computed inv(A) if possible

    function GainNode(A::Union{Array{Float64},Float64}=1.0; id=generateNodeId(GainNode))
        # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
        A = (typeof(A)==Float64) ? fill!(Array(Float64,1,1),A) : ensureMatrix(deepcopy(A))
        self = new(A, id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(current_graph, self)

        for (iface_index, iface_handle) in enumerate([:in, :out, :gain])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        # Try to precompute inv(A)
        try
            self.A_inv = inv(self.A)
        catch
            warn("The specified multiplier for $(typeof(self)) $(self.id) is not invertible. This might cause problems. Double check that this is what you want.")
        end

        return self
    end

end

isDeterministic(::GainNode) = true

############################################
# GaussianDistribution methods
############################################

# Backward Gaussian to IN
function sumProduct!(   node::GainNode,
                        outbound_interface_index::Int,
                        ::Void,
                        msg_out::Message{GaussianDistribution})
    (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
    dist_2 = ensureParameters!(msg_out.payload, (:xi, :W))
    dist_1.xi = node.A[1,1] * dist_2.xi
    dist_1.W = (node.A[1,1])^2 * dist_2.W
    dist_1.m = dist_1.V = NaN

    return (:gain_gaussian_backward,
            node.interfaces[outbound_interface_index].message)
end

# Forward Gaussian to OUT
function sumProduct!(   node::GainNode,
                        outbound_interface_index::Int,
                        msg_in::Message{GaussianDistribution},
                        ::Void)
    (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
    dist_1 = ensureParameters!(msg_in.payload, (:m, :V))
    dist_2.m = node.A[1,1] * dist_1.m
    dist_2.V = (node.A[1,1])^2 * dist_1.V
    dist_2.xi = dist_2.W = NaN

    return (:gain_gaussian_forward,
            node.interfaces[outbound_interface_index].message)
end

#Forward Gaussian to OUT if gain is present on the edge
function sumProduct!(   node::GainNode,
                        outbound_interface_index::Int,
                        msg_in::Message{GaussianDistribution},
                        msg_gain::Message{DeltaDistribution},
                        ::Void)
  (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")

  dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
  dist_1 = ensureParameters!(msg_in.payload, (:m, :V))
  dist_2.m = gain.m * dist_1.m
  dist_2.V = (gain.m)^2 * dist_1.V
  dist_2.xi = dist_2.W = NaN

  return (:gain_gaussian_forward,
          node.interfaces[outbound_interface_index].message)
end

#Backward Gaussian to IN if gain is present on the edge
function sumProduct!(   node::GainNode,
                        outbound_interface_index::Int,
                        msg_in::Message{GaussianDistribution},
                        msg_gain::Message{DeltaDistribution},
                        ::Void)
  (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
  dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
  dist_2 = ensureParameters!(msg_out.payload, (:xi, :W))
  dist_1.xi = gain.m * dist_2.xi
  dist_1.W = (gain.m)^2 * dist_2.W
  dist_1.m = dist_1.V = NaN

  return (:gain_gaussian_backward,
          node.interfaces[outbound_interface_index].message)
end

############################################
# DeltaDistribution methods
############################################

# Backward DeltaDistribution to IN
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_out::Message{DeltaDistribution{Float64}})
    (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload
    dist_1.m = node.A_inv[1,1] * msg_out.payload.m

    return (:gain_delta_backward,
            node.interfaces[outbound_interface_index].message)
end

# Forward DeltaDistribution to OUT
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     msg_in::Message{DeltaDistribution{Float64}},
                     ::Void)
    (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload
    dist_2.m = node.A[1,1] * msg_in.payload.m

    return (:gain_delta_forward,
            node.interfaces[outbound_interface_index].message)
end

############################################
# MvGaussianDistribution methods
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
function fixedGainGaussianBackwardRule!(dist_result::MvGaussianDistribution, dist_2::MvGaussianDistribution, A::Any, A_inv::Any=nothing)
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
        ensureParameters!(dist_2, (:xi, :W))
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardFixedGainWRule(A, dist_2.W)
        dist_result.xi = backwardFixedGainXiRule(A, dist_2.xi)
    end

    return dist_result
end

function sumProduct!(node::GainNode,
                            outbound_interface_index::Int,
                            ::Void,
                            msg_out::Message{MvGaussianDistribution})
    isProper(msg_out.payload) || error("Improper input distributions are not supported")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], MvGaussianDistribution).payload

    if outbound_interface_index == 1
        fixedGainGaussianBackwardRule!(dist_1, msg_out.payload, node.A, isdefined(node, :A_inv) ? node.A_inv : nothing)
    else
        error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    end

    return (:fixed_gain_gaussian_backward,
            node.interfaces[outbound_interface_index].message)
end

# Forward Gaussian to OUT
function fixedGainGaussianForwardRule!(dist_result::MvGaussianDistribution, dist_1::MvGaussianDistribution, A::Any, A_inv::Any=nothing)
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
        ensureParameters!(dist_1, (:m, :V))
        dist_result.m = forwardFixedGainMRule(A, dist_1.m)
        dist_result.V = forwardFixedGainVRule(A, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    end

    return dist_result
end

function sumProduct!(node::GainNode,
                            outbound_interface_index::Int,
                            msg_in::Message{MvGaussianDistribution},
                            ::Void)
    isProper(msg_in.payload) || error("Improper input distributions are not supported")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], MvGaussianDistribution).payload

    if outbound_interface_index == 2
        # Forward message
        dist_1 = msg_in.payload

        fixedGainGaussianForwardRule!(dist_2, msg_in.payload, node.A, isdefined(node, :A_inv) ? node.A_inv : nothing)
    else
        error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    end

    return (:fixed_gain_gaussian_forward,
            node.interfaces[outbound_interface_index].message)
end


############################################
# MvDeltaDistribution methods
############################################

# Backward MvDeltaDistribution to IN
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_out::Message{MvDeltaDistribution{Float64}})
    (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], MvDeltaDistribution{Float64}).payload
    dist_1.m = node.A_inv * msg_out.payload.m

    return (:fixed_gain_delta_backward,
            node.interfaces[outbound_interface_index].message)
end

# Forward MvDeltaDistribution to OUT
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     msg_in::Message{MvDeltaDistribution{Float64}},
                     ::Void)
    (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], MvDeltaDistribution{Float64}).payload
    dist_2.m = node.A * msg_in.payload.m

    return (:fixed_gain_delta_forward,
            node.interfaces[outbound_interface_index].message)
end
