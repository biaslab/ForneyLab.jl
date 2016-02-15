export GainEqualityNode

"""
Description:

    Gain-equality node: A⁻¹*out = in1 = in2
    Combines the node functions of the GainNode
    and the EqualityNode for computational efficiency.

         _________
     in1 |       | in2
    -----|->[=]<-|-----
         |   |   |
         |   v   |
         |  [A]  |
         |___|___|
             | out
             v

    f(in1,in2,out) = δ(A*in1 - out)⋅δ(A*in2 - out)

Construction:

    GainEqualityNode([1.0], id=:my_node)
"""
type GainEqualityNode <: Node
    A::Array{Float64}
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible

    function GainEqualityNode(A::Union{Array{Float64},Float64}=1.0; id=generateNodeId(GainEqualityNode))
        self = new(ensureMatrix(deepcopy(A)), id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:in1, :in2, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        # Precompute inverse of A
        try self.A_inv = inv(self.A) end

        return self
    end
end

isDeterministic(::GainEqualityNode) = true


############################################
# GaussianDistribution methods
############################################

"""
GainEqualityNode:

        _________
      N |       | N
    ----|->[=]<-|----
        |   |   |
        |   v   |
        |  [A]  |
        |___|___|
          | | N
          v v

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainEqualityNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Any)

    dist_temp = GaussianDistribution()
    equalityRule!(dist_temp, msg_in1.payload, msg_in2.payload)
    dist_temp = ensureParameters!(dist_temp, (:m, :V))

    outbound_dist.m = node.A[1,1] * dist_temp.m
    outbound_dist.V = (node.A[1,1])^2 * dist_temp.V
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GainEqualityNode:

        _________
      N |       | N
    ----|->[=]<-|----
        |   |   | -->
        |   v   |
        |  [A]  |
        |___|___|
            | N
            v

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainEqualityNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Any,
                            msg_out::Message{GaussianDistribution})

    return gainEqualityBackwardRule!(outbound_dist, msg_in1.payload, msg_out.payload, node.A)
end

"""
GainEqualityNode:

        _________
      N |       | N
    ----|->[=]<-|----
    <-- |   |   |
        |   v   |
        |  [A]  |
        |___|___|
            | N
            v

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainEqualityNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_in1::Void,
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Message{GaussianDistribution})
    # Backward message (towards in1)
    return gainEqualityBackwardRule!(outbound_dist, msg_in2.payload, msg_out.payload, node.A)
end

function gainEqualityBackwardRule!(dist_result::GaussianDistribution, dist_in::GaussianDistribution, dist_out::GaussianDistribution, A::Any)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    # Backward message (towards in1 or in2)
    ensureParameters!(dist_out, (:xi, :W))
    ensureParameters!(dist_in, (:xi, :W))

    dist_result.m = NaN
    dist_result.V = NaN
    dist_result.W = dist_in.W + A[1,1]^2 * dist_out.W
    dist_result.xi = dist_in.xi + A[1,1] * dist_out.xi

    return dist_result
end


############################################
# MvGaussianDistribution methods
############################################

function sumProductRule!{T<:MvGaussianDistribution}(node::GainEqualityNode,
                                                    outbound_interface_index::Type{Val{3}},
                                                    outbound_dist::T,
                                                    msg_in1::Message{T},
                                                    msg_in2::Message{T},
                                                    msg_out::Any)

    dist_temp = vague(MvGaussianDistribution{size(node.A, 2)})
    equalityRule!(dist_temp, msg_in1.payload, msg_in2.payload)
    return gainForwardRule!(outbound_dist, dist_temp, node.A, (isdefined(node, :A_inv)) ? node.A_inv : nothing)
end

function sumProductRule!{T<:MvGaussianDistribution}(node::GainEqualityNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::T,
                                                    msg_in1::Message{T},
                                                    msg_in2::Any,
                                                    msg_out::Message{T})

    return gainEqualityBackwardRule!(outbound_dist, msg_in1.payload, msg_out.payload, node.A)
end

function sumProductRule!{T<:MvGaussianDistribution}(node::GainEqualityNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::T,
                                                    msg_in1::Void,
                                                    msg_in2::Message{T},
                                                    msg_out::Message{T})

    return gainEqualityBackwardRule!(outbound_dist, msg_in2.payload, msg_out.payload, node.A)
end

function gainEqualityBackwardRule!( dist_result::MvGaussianDistribution,
                                    dist_in::MvGaussianDistribution,
                                    dist_out::MvGaussianDistribution,
                                    A::Any)

    # Select parameterization
    # Order is from least to most computationally intensive
    if isValid(dist_out.xi) && isValid(dist_out.W) && isValid(dist_in.xi) && isValid(dist_in.W)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardGainEqualityWRule(A, dist_in.W, dist_out.W)
        dist_result.xi = backwardGainEqualityXiRule(A, dist_in.xi, dist_out.xi)
    elseif isValid(dist_out.m) && isValid(dist_out.V) && isValid(dist_in.m) && isValid(dist_in.V)
        dist_result.m = backwardGainEqualityMRule(A, dist_in.m, dist_in.V, dist_out.m, dist_out.V)
        dist_result.V = backwardGainEqualityVRule(A, dist_in.V, dist_out.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_out.m) && isValid(dist_out.W) && isValid(dist_in.m) && isValid(dist_in.W)
        dist_result.m = backwardGainEqualityMRule(A, dist_in.m, inv(dist_in.W), dist_out.m, inv(dist_out.W))
        invalidate!(dist_result.V)
        dist_result.W = backwardGainEqualityWRule(A, dist_in.W, dist_out.W)
        invalidate!(dist_result.xi)
    else
        # Fallback: convert inbound messages to (xi,W) parametrization and then use efficient rules
        ensureParameters!(dist_in, (:xi, :W))
        ensureParameters!(dist_out, (:xi, :W))
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardGainEqualityWRule(A, dist_in.W, dist_out.W)
        dist_result.xi = backwardGainEqualityXiRule(A, dist_in.xi, dist_out.xi)
    end

    return dist_result
end

# Rule set for backward propagation ({in2,out}-->in1 or {in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y
backwardGainEqualityVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x - V_x * A' * inv(V_y + A * V_x * A') * A * V_x
backwardGainEqualityMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, V_x::Array{T, 2}, m_y::Array{T, 1}, V_y::Array{T, 2}) = m_x + V_x * A' * inv(V_y + A * V_x * A') * (m_y - A * m_x)
