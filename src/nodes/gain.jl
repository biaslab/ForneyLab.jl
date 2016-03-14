export GainNode

"""
Description:

    Multiplication:
        out = gain * in

        gain
         |
    in   V   out
    ----->[A]----->

    f(in,out,gain) =  δ(out - gain*in), where gain is either provided upon construction of the node and is a fixed value or is supplied via gain interface.

Interfaces:

    1 i[:in], 2 i[:out], 3 i[:gain] (optional)

Construction:

    GainNode(gain=[1.0], id=:my_node)
    or
    GainNode(id=:my_node)
"""
type GainNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    gain::AbstractMatrix{Float64}
    gain_inv::AbstractMatrix{Float64} # holds pre-computed inv(gain) if possible

    function GainNode(;gain::Union{AbstractArray{Float64}, Float64, Void}=nothing, id=generateNodeId(GainNode))
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

        addNode!(currentGraph(), self)

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

"""
GainNode:

      N        N
    --->[gain]--->
    <--

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_in::Any,
                            msg_out::Message{GaussianDistribution})

    dist_out = ensureParameters!(msg_out.payload, (:xi, :W))

    outbound_dist.m = NaN
    outbound_dist.V = NaN
    outbound_dist.xi = node.gain[1,1] * dist_out.xi
    outbound_dist.W = (node.gain[1,1])^2 * dist_out.W

    return outbound_dist
end

"""
GainNode:

      N        N
    --->[gain]--->
               -->

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_in::Message{GaussianDistribution},
                            msg_out::Any)

    dist_in = ensureParameters!(msg_in.payload, (:m, :V))

    outbound_dist.m = node.gain[1,1] * dist_in.m
    outbound_dist.V = (node.gain[1,1])^2 * dist_in.V
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GainNode:

          | δ
      N   v    N
    --->[gain]--->
               -->

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_in::Message{GaussianDistribution},
                            msg_out::Any,
                            msg_gain::Message{DeltaDistribution{Float64}})

    dist_in = ensureParameters!(msg_in.payload, (:m, :V))
    gain_dist = msg_gain.payload

    outbound_dist.m = gain_dist.m * dist_in.m
    outbound_dist.V = (gain_dist.m)^2 * dist_in.V
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GainNode:

          | δ
      N   v    N
    --->[gain]--->
    <--

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_in::Any,
                            msg_out::Message{GaussianDistribution},
                            msg_gain::Message{DeltaDistribution{Float64}})

    dist_out = ensureParameters!(msg_out.payload, (:xi, :W))
    gain_dist = msg_gain.payload

    outbound_dist.m = NaN
    outbound_dist.V = NaN
    outbound_dist.xi = gain_dist.m * dist_out.xi
    outbound_dist.W = (gain_dist.m)^2 * dist_out.W

    return outbound_dist
end


############################################
# DeltaDistribution methods
############################################

"""
GainNode:

      δ        δ
    --->[gain]--->
    <--
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::DeltaDistribution{Float64},
                            msg_in::Any,
                            msg_out::Message{DeltaDistribution{Float64}})

    outbound_dist.m = node.gain_inv[1,1] * msg_out.payload.m
    return outbound_dist
end

"""
GainNode:

          | δ
      δ   v    δ
    --->[gain]--->
    <--
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::DeltaDistribution{Float64},
                            msg_in::Any,
                            msg_out::Message{DeltaDistribution{Float64}},
                            msg_gain::Message{DeltaDistribution{Float64}})

    outbound_dist.m = msg_out.payload.m / msg_gain.payload.m
    return outbound_dist
end

"""
GainNode:

      δ        δ
    --->[gain]--->
               -->
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::DeltaDistribution{Float64},
                            msg_in::Message{DeltaDistribution{Float64}},
                            msg_out::Any)

    outbound_dist.m = node.gain[1,1] * msg_in.payload.m
    return outbound_dist
end

"""
GainNode:

          | δ
      δ   v    δ
    --->[gain]--->
               -->
"""
function sumProductRule!(   node::GainNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::DeltaDistribution{Float64},
                            msg_in::Message{DeltaDistribution{Float64}},
                            msg_out::Any,
                            msg_gain::Message{DeltaDistribution{Float64}})

    outbound_dist.m = msg_gain.payload.m * msg_in.payload.m
    return outbound_dist
end


############################################
# MvGaussianDistribution methods
############################################

function sumProductRule!{T<:MvGaussianDistribution}(node::GainNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::T,
                                                    msg_in::Any,
                                                    msg_out::Message{T})

    return gainBackwardRule!(outbound_dist, msg_out.payload, node.gain, isdefined(node, :gain_inv) ? node.gain_inv : nothing)
end

function sumProductRule!{T<:MvGaussianDistribution}(node::GainNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::T,
                                                    msg_in::Message{T},
                                                    msg_out::Any)

    return gainForwardRule!(outbound_dist, msg_in.payload, node.gain, isdefined(node, :gain_inv) ? node.gain_inv : nothing)
end

function sumProductRule!{T<:MvGaussianDistribution}(node::GainNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::T,
                                                    msg_in::Any,
                                                    msg_out::Message{T},
                                                    msg_gain::Message{DeltaDistribution{Matrix{Float64}}})

    return gainBackwardRule!(outbound_dist, msg_out.payload, msg_gain.payload.m)
end

function sumProductRule!{T<:MvGaussianDistribution}(node::GainNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::T,
                                                    msg_in::Message{T},
                                                    msg_out::Any,
                                                    msg_gain::Message{DeltaDistribution{Matrix{Float64}}})

    return gainForwardRule!(outbound_dist, msg_in.payload, msg_gain.payload.m)
end

function gainBackwardRule!(dist_result::MvGaussianDistribution, dist_2::MvGaussianDistribution, A::Any, A_inv::Any=nothing)
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

function gainForwardRule!(dist_result::MvGaussianDistribution, dist_1::MvGaussianDistribution, A::Any, gain_inv::Any=nothing)
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

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainMRule{T<:Number}(A_inv::AbstractMatrix{T}, m::Vector{T}) = A_inv * m
backwardGainMRule{T<:Number}(A::AbstractMatrix{T}, m::Vector{T}, W::AbstractMatrix{T}) = pinv(A' * W * A) * A' * W * m
backwardGainVRule{T<:Number}(A_inv::AbstractMatrix{T}, V::AbstractMatrix{T}) = A_inv * V * A_inv'
backwardGainWRule{T<:Number}(A::AbstractMatrix{T}, W::AbstractMatrix{T}) = A' * W * A
backwardGainXiRule{T<:Number}(A::AbstractMatrix{T}, xi::Vector{T}) = A' * xi

# Rule set for forward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainMRule{T<:Number}(A::AbstractMatrix{T}, m::Vector{T}) = A * m
forwardGainVRule{T<:Number}(A::AbstractMatrix{T}, V::AbstractMatrix{T}) = A * V * A'
forwardGainWRule{T<:Number}(A_inv::AbstractMatrix{T}, W::AbstractMatrix{T}) = A_inv' * W * A_inv
forwardGainXiRule{T<:Number}(A_inv::AbstractMatrix{T}, xi::Vector{T}) = A_inv' * xi
forwardGainXiRule{T<:Number}(A::AbstractMatrix{T}, xi::Vector{T}, V::AbstractMatrix{T}) = pinv(A * V * A') * A * V * xi # Combination of xi and V


############################################
# MvDeltaDistribution methods
############################################

function sumProductRule!{T<:MvDeltaDistribution{Float64}}(  node::GainNode,
                                                            outbound_interface_index::Type{Val{1}},
                                                            outbound_dist::T,
                                                            msg_in::Any,
                                                            msg_out::Message{T})

    outbound_dist.m = node.gain_inv * msg_out.payload.m
    return outbound_dist
end

function sumProductRule!{T<:MvDeltaDistribution{Float64}}(  node::GainNode,
                                                            outbound_interface_index::Type{Val{1}},
                                                            outbound_dist::T,
                                                            msg_in::Any,
                                                            msg_out::Message{T},
                                                            msg_gain::Message{DeltaDistribution{Matrix{Float64}}})

    outbound_dist.m = inv(msg_gain.payload.m) * msg_out.payload.m
    return outbound_dist
end

function sumProductRule!{T<:MvDeltaDistribution{Float64}}(  node::GainNode,
                                                            outbound_interface_index::Type{Val{2}},
                                                            outbound_dist::T,
                                                            msg_in::Message{T},
                                                            msg_out::Any)

    outbound_dist.m = node.gain * msg_in.payload.m
    return outbound_dist
end

function sumProductRule!{T<:MvDeltaDistribution{Float64}}(  node::GainNode,
                                                            outbound_interface_index::Type{Val{2}},
                                                            outbound_dist::T,
                                                            msg_in::Message{T},
                                                            msg_out::Any,
                                                            msg_gain::Message{DeltaDistribution{Matrix{Float64}}})

    outbound_dist.m = msg_gain.payload.m * msg_in.payload.m
    return outbound_dist
end
