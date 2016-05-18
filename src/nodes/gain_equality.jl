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
    gain::AbstractMatrix{Float64}
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    gain_inv::AbstractMatrix{Float64} # holds pre-computed inv(gain) if possible

    function GainEqualityNode(gain::Union{AbstractArray{Float64}, Float64}=1.0; id=generateNodeId(GainEqualityNode))
        self = new(ensureMatrix(deepcopy(gain)), id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:in1, :in2, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        # Precompute inverse of gain
        self.gain_inv = pinv(self.gain)

        return self
    end
end

isDeterministic(::GainEqualityNode) = true

outboundParameterValue(node::GainEqualityNode, ::Type{Val{:dims_n}}, args...) = size(node.gain, 1)

outboundParameterValue(node::GainEqualityNode, ::Type{Val{:dims_m}}, args...) = size(node.gain, 2)

############################################
# Gaussian methods
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
                            outbound_dist::Gaussian,
                            msg_in1::Message{Gaussian},
                            msg_in2::Message{Gaussian},
                            msg_out::Any)

    dist_temp = ensureParameters!(msg_in1.payload * msg_in2.payload, (:m, :V))

    outbound_dist.m = node.gain[1,1] * dist_temp.m
    outbound_dist.V = (node.gain[1,1])^2 * dist_temp.V
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
                            outbound_dist::Gaussian,
                            msg_in1::Message{Gaussian},
                            msg_in2::Any,
                            msg_out::Message{Gaussian})

    return gainEqualityBackwardRule!(outbound_dist, msg_in1.payload, msg_out.payload, node.gain)
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
                            outbound_dist::Gaussian,
                            msg_in1::Void,
                            msg_in2::Message{Gaussian},
                            msg_out::Message{Gaussian})
    # Backward message (towards in1)
    return gainEqualityBackwardRule!(outbound_dist, msg_in2.payload, msg_out.payload, node.gain)
end

function gainEqualityBackwardRule!(dist_result::Gaussian, dist_in::Gaussian, dist_out::Gaussian, A::Any)
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
# MvGaussian methods
############################################

function sumProductRule!{dims_n, dims_m}(   node::GainEqualityNode,
                                            outbound_interface_index::Type{Val{3}},
                                            outbound_dist::MvGaussian{dims_n},
                                            msg_in1::Message{MvGaussian{dims_m}},
                                            msg_in2::Message{MvGaussian{dims_m}},
                                            msg_out::Any)

    dist_temp = msg_in1.payload * msg_in2.payload
    return gainForwardRule!(outbound_dist, dist_temp, node.gain, (isdefined(node, :A_inv)) ? node.gain_inv : nothing)
end

function sumProductRule!{dims_n, dims_m}(   node::GainEqualityNode,
                                            outbound_interface_index::Type{Val{2}},
                                            outbound_dist::MvGaussian{dims_m},
                                            msg_in1::Message{MvGaussian{dims_m}},
                                            msg_in2::Any,
                                            msg_out::Message{MvGaussian{dims_n}})

    return gainEqualityBackwardRule!(outbound_dist, msg_in1.payload, msg_out.payload, node.gain)
end

function sumProductRule!{dims_n, dims_m}(   node::GainEqualityNode,
                                            outbound_interface_index::Type{Val{1}},
                                            outbound_dist::MvGaussian{dims_m},
                                            msg_in1::Any,
                                            msg_in2::Message{MvGaussian{dims_m}},
                                            msg_out::Message{MvGaussian{dims_n}})

    return gainEqualityBackwardRule!(outbound_dist, msg_in2.payload, msg_out.payload, node.gain)
end

function gainEqualityBackwardRule!( dist_result::MvGaussian,
                                    dist_in::MvGaussian,
                                    dist_out::MvGaussian,
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
        dist_result.m = backwardGainEqualityMRule(A, dist_in.m, cholinv(dist_in.W), dist_out.m, cholinv(dist_out.W))
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
backwardGainEqualityWRule{T<:Number}(A::AbstractMatrix{T}, W_x::AbstractMatrix{T}, W_y::AbstractMatrix{T}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::AbstractMatrix{T}, xi_x::Vector{T}, xi_y::Vector{T}) = xi_x + A' * xi_y
backwardGainEqualityVRule{T<:Number}(A::AbstractMatrix{T}, V_x::AbstractMatrix{T}, V_y::AbstractMatrix{T}) = V_x - V_x * A' * cholinv(V_y + A * V_x * A') * A * V_x
backwardGainEqualityMRule{T<:Number}(A::AbstractMatrix{T}, m_x::Vector{T}, V_x::AbstractMatrix{T}, m_y::Vector{T}, V_y::AbstractMatrix{T}) = m_x + V_x * A' * cholinv(V_y + A * V_x * A') * (m_y - A * m_x)
