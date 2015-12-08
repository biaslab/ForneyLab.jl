############################################
# GainNode
############################################
# Description:
#   Multiplication:
#         out = gain * in
#
#         gain
#          |
#     in   V   out
#   ----->[A]----->
#
#
#
#     f(in,out,gain) =  δ(out - gain*in), where gain is either provided upon construction of the node and is a fixed value or is supplied via gain interface.
#
# Interfaces:
#   1 i[:in], 2 i[:out], 3 i[:gain] (optional)
#
# Construction:
#    GainNode(gain=[1.0], id=:my_node)
#     or
#    GainNode(id=:my_node)
#
############################################

export GainNode

type GainNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    gain::Matrix{Float64}
    gain_inv::Matrix{Float64} # holds pre-computed inv(gain) if possible

    function GainNode(;gain::Union{Array{Float64},Float64, Void}=nothing, id=generateNodeId(GainNode))
        if gain != nothing
            # Fixed gain; no gain interface.
            # Deepcopy gain to avoid an unexpected change of the input argument gain. Ensure that gain is a matrix.
            gain = (typeof(gain)==Float64) ? fill!(Array(Float64,1,1),gain) : ensureMatrix(deepcopy(gain))
            self = new(id, Array(Interface, 2), Dict{Symbol,Interface}(), gain)
            # Try to precompute inv(gain)
            try
                self.gain_inv = inv(self.gain)
            catch
                warn("The specified gain for $(typeof(self)) $(self.id) is not invertible. This might cause problems. Double check that this is what you want.")
            end
        else
            self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        end

        addNode!(current_graph, self)

        for (iface_index, iface_handle) in enumerate([:in, :out, :gain][1:length(self.interfaces)])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
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
    dist_1.xi = node.gain[1,1] * dist_2.xi
    dist_1.W = (node.gain[1,1])^2 * dist_2.W
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
    dist_2.m = node.gain[1,1] * dist_1.m
    dist_2.V = (node.gain[1,1])^2 * dist_1.V
    dist_2.xi = dist_2.W = NaN

    return (:gain_gaussian_forward,
            node.interfaces[outbound_interface_index].message)
end

# Forward Gaussian to OUT if gain is present on the edge
function sumProduct!(   node::GainNode,
                        outbound_interface_index::Int,
                        msg_in::Message{GaussianDistribution},
                        ::Void,
                        msg_gain::Message{DeltaDistribution{Float64}})
  (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")

  dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
  dist_1 = ensureParameters!(msg_in.payload, (:m, :V))
  gain_dist = msg_gain.payload

  dist_2.m = gain_dist.m * dist_1.m
  dist_2.V = (gain_dist.m)^2 * dist_1.V
  dist_2.xi = dist_2.W = NaN

  return (:gain_gaussian_forward,
          node.interfaces[outbound_interface_index].message)
end

# Backward Gaussian to IN if gain is present on the edge
function sumProduct!(   node::GainNode,
                        outbound_interface_index::Int,
                        ::Void,
                        msg_out::Message{GaussianDistribution},
                        msg_gain::Message{DeltaDistribution{Float64}})

  (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
  dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
  dist_2 = ensureParameters!(msg_out.payload, (:xi, :W))
  gain_dist = msg_gain.payload

  dist_1.xi = gain_dist.m * dist_2.xi
  dist_1.W = (gain_dist.m)^2 * dist_2.W
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
    dist_1.m = node.gain_inv[1,1] * msg_out.payload.m

    return (:gain_delta_backward,
            node.interfaces[outbound_interface_index].message)
end

function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_out::Message{DeltaDistribution{Float64}},
                     msg_gain::Message{DeltaDistribution{Float64}})
    (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload
    dist_1.m = msg_out.payload.m / msg_gain.payload.m

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
    dist_2.m = node.gain[1,1] * msg_in.payload.m

    return (:gain_delta_forward,
            node.interfaces[outbound_interface_index].message)
end

function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     msg_in::Message{DeltaDistribution{Float64}},
                     ::Void,
                     msg_gain::Message{DeltaDistribution{Float64}})
    (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload
    dist_2.m = msg_gain.payload.m * msg_in.payload.m

    return (:gain_delta_forward,
            node.interfaces[outbound_interface_index].message)
end

############################################
# MvGaussianDistribution methods
############################################

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainMRule{T<:Number}(A_inv::Array{T, 2}, m::Array{T, 1}) = A_inv * m
backwardGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}, W::Array{T, 2}) = pinv(A' * W * A) * A' * W * m
backwardGainVRule{T<:Number}(A_inv::Array{T, 2}, V::Array{T, 2}) = A_inv * V * A_inv'
backwardGainWRule{T<:Number}(A::Array{T, 2}, W::Array{T, 2}) = A' * W * A
backwardGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}) = A' * xi

# Rule set for forward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainMRule{T<:Number}(A::Array{T, 2}, m::Array{T, 1}) = A * m
forwardGainVRule{T<:Number}(A::Array{T, 2}, V::Array{T, 2}) = A * V * A'
forwardGainWRule{T<:Number}(A_inv::Array{T, 2}, W::Array{T, 2}) = A_inv' * W * A_inv
forwardGainXiRule{T<:Number}(A_inv::Array{T, 2}, xi::Array{T, 1}) = A_inv' * xi
forwardGainXiRule{T<:Number}(A::Array{T, 2}, xi::Array{T, 1}, V::Array{T, 2}) = pinv(A * V * A') * A * V * xi # Combination of xi and V

# Backward Gaussian to IN
function gainGaussianBackwardRule!(dist_result::MvGaussianDistribution, dist_2::MvGaussianDistribution, A::Any, A_inv::Any=nothing)
    # Calculations for a gaussian message type; Korl (2005), table 4.1
    # dist_result = inv(A) * dist_2

    if isValid(dist_2.xi) && isValid(dist_2.W)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardGainWRule(A, dist_2.W)
        dist_result.xi = backwardGainXiRule(A, dist_2.xi)
    elseif isValid(dist_2.m) && isValid(dist_2.W) && isRoundedPosDef(A) && isRoundedPosDef(dist_2.W) && (A_inv != nothing)
        dist_result.m = backwardGainMRule(A_inv, dist_2.m) # Short version of the rule
        invalidate!(dist_result.V)
        dist_result.W = backwardGainWRule(A, dist_2.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_2.m) && isValid(dist_2.W)
        dist_result.m = backwardGainMRule(A, dist_2.m, dist_2.W) # Long version of the rule
        invalidate!(dist_result.V)
        dist_result.W = backwardGainWRule(A, dist_2.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_2.m) && isValid(dist_2.V) && (A_inv != nothing) && (dist_2.V==zeros(size(dist_2.V)) || isRoundedPosDef(dist_2.V)) # V positive definite <=> W (its inverse) is also positive definite. Border-case is an all-zero V, in which case the outbound message also has zero variance.
        dist_result.m = backwardGainMRule(A_inv, dist_2.m) # Short version of the rule, only valid if A is positive definite and (W is positive definite or V is 0)
        dist_result.V = backwardGainVRule(A_inv, dist_2.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    else
        # Fallback: convert inbound message to (xi,W) parametrization and then use efficient rules
        ensureParameters!(dist_2, (:xi, :W))
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardGainWRule(A, dist_2.W)
        dist_result.xi = backwardGainXiRule(A, dist_2.xi)
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
        gainGaussianBackwardRule!(dist_1, msg_out.payload, node.gain, isdefined(node, :gain_inv) ? node.gain_inv : nothing)
    else
        error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    end

    return (:gain_gaussian_backward,
            node.interfaces[outbound_interface_index].message)
end

# Forward Gaussian to OUT
function gainGaussianForwardRule!(dist_result::MvGaussianDistribution, dist_1::MvGaussianDistribution, A::Any, gain_inv::Any=nothing)
    # Calculations for a gaussian message type; Korl (2005), table 4.1
    # dist_result = A * dist_1

    if isValid(dist_1.m) && isValid(dist_1.V)
        dist_result.m = forwardGainMRule(A, dist_1.m)
        dist_result.V = forwardGainVRule(A, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.m) && isValid(dist_1.W) && (gain_inv != nothing)
        dist_result.m = forwardGainMRule(A, dist_1.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardGainWRule(gain_inv, dist_1.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.xi) && isValid(dist_1.W) && isRoundedPosDef(A) && isRoundedPosDef(dist_1.W) && (gain_inv != nothing) # V should be positive definite, meaning V is invertible and its inverse W is also positive definite.
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardGainWRule(gain_inv, dist_1.W)
        dist_result.xi = forwardGainXiRule(gain_inv, dist_1.xi) # Short version of the rule, only valid if A and V are positive definite
    elseif isValid(dist_1.xi) && isValid(dist_1.V) && isValid(dist_1.W) && (gain_inv != nothing)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardGainWRule(gain_inv, dist_1.W)
        dist_result.xi = forwardGainXiRule(gain_inv, dist_1.xi, dist_1.V) # Long version of the rule
    else
        # Fallback: convert inbound message to (m,V) parametrization and then use efficient rules
        ensureParameters!(dist_1, (:m, :V))
        dist_result.m = forwardGainMRule(A, dist_1.m)
        dist_result.V = forwardGainVRule(A, dist_1.V)
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
        gainGaussianForwardRule!(dist_2, msg_in.payload, node.gain, isdefined(node, :gain_inv) ? node.gain_inv : nothing)
    else
        error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    end

    return (:gain_gaussian_forward,
            node.interfaces[outbound_interface_index].message)
end

# Backward Gaussian to IN if gain is present on the edge
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_out::Message{MvGaussianDistribution},
                     msg_gain::Message{DeltaDistribution{Matrix{Float64}}})

    isProper(msg_out.payload) || error("Improper input distributions are not supported")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], MvGaussianDistribution).payload

    if outbound_interface_index == 1
        # We don't pass the inverse gain because we can't cache it if the gain comes in as a message.
        gainGaussianBackwardRule!(dist_1, msg_out.payload, msg_gain.payload.m)
     else
        error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
     end

     return (:gain_gaussian_forward,
             node.interfaces[outbound_interface_index].message)
end

# Forward Gaussian to OUT if gain is present on the edge
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     msg_in::Message{MvGaussianDistribution},
                     ::Void,
                     msg_gain::Message{DeltaDistribution{Matrix{Float64}}})

     isProper(msg_in.payload) || error("Improper input distributions are not supported")
     dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], MvGaussianDistribution).payload

     if outbound_interface_index == 2
         # We don't pass the inverse gain because we can't cache it if the gain comes in as a message.
         gainGaussianForwardRule!(dist_2, msg_in.payload, msg_gain.payload.m)
     else
         error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
     end

     return (:gain_gaussian_forward,
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
    dist_1.m = node.gain_inv * msg_out.payload.m

    return (:gain_delta_backward,
            node.interfaces[outbound_interface_index].message)
end

function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_out::Message{MvDeltaDistribution{Float64}},
                     msg_gain::Message{DeltaDistribution{Matrix{Float64}}})
    (outbound_interface_index == 1) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_1 = ensureMessage!(node.interfaces[outbound_interface_index], MvDeltaDistribution{Float64}).payload
    dist_1.m = inv(msg_gain.payload.m) * msg_out.payload.m

    return (:gain_delta_backward,
            node.interfaces[outbound_interface_index].message)
end

# Forward MvDeltaDistribution to OUT
function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     msg_in::Message{MvDeltaDistribution{Float64}},
                     ::Void)
    (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], MvDeltaDistribution{Float64}).payload
    dist_2.m = node.gain * msg_in.payload.m

    return (:gain_delta_forward,
            node.interfaces[outbound_interface_index].message)
end

function sumProduct!(node::GainNode,
                     outbound_interface_index::Int,
                     msg_in::Message{MvDeltaDistribution{Float64}},
                     ::Void,
                     msg_gain::Message{DeltaDistribution{Matrix{Float64}}})
    (outbound_interface_index == 2) || error("Invalid interface id $(outbound_interface_index) for calculating message on $(typeof(node)) $(node.id)")
    dist_2 = ensureMessage!(node.interfaces[outbound_interface_index], MvDeltaDistribution{Float64}).payload
    dist_2.m = msg_gain.payload.m * msg_in.payload.m

    return (:gain_delta_forward,
            node.interfaces[outbound_interface_index].message)
end
