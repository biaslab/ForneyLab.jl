export GainAdditionNode

"""
Description:

    Gain-addition node: out = A*in1 + in2
    Combines the node functions of the GainNode
    and the AdditionNode for computational efficiency.

             | in1
             |
         ____|____
         |   v   |
         |  [A]  |
         |   |   |
     in2 |   v   | out
    -----|->[+]--|---->
         |_______|

    f(in1,in2,out) = δ(out - A*in1 - in2)

Interfaces:

    1 i[:in1], 2 i[:in2], 3 i[:out]

Construction:

    GainAdditionNode([1.0], id=:my_node)
"""
type GainAdditionNode <: Node
    gain::AbstractMatrix{Float64}
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    gain_inv::AbstractMatrix{Float64} # holds pre-computed inv(gain) if possible

    function GainAdditionNode(gain::Union{AbstractArray{Float64}, Float64}=1.0; id=generateNodeId(GainAdditionNode))
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

isDeterministic(::GainAdditionNode) = true

outboundParameterValue(node::GainAdditionNode, ::Type{Val{:dims_n}}, args...) = size(node.gain, 1)

outboundParameterValue(node::GainAdditionNode, ::Type{Val{:dims_m}}, args...) = size(node.gain, 2)

############################################
# Gaussian methods
############################################

"""
GainAdditionNode:

            | N
        ____|____
        |   v   |
        |  [A]  |
        |   |   |
      N |   v   | N
    ----|->[+]--|--->
        |_______| -->

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainAdditionNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            msg_in1::Message{Gaussian},
                            msg_in2::Message{Gaussian},
                            msg_out::Any)

    dist_1 = ensureParameters!(msg_in1.payload, (:m, :V))
    dist_2 = ensureParameters!(msg_in2.payload, (:m, :V))

    outbound_dist.m = dist_2.m + node.gain[1,1]*dist_1.m
    outbound_dist.V = dist_2.V + node.gain[1,1]^2 * dist_1.V
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GainAdditionNode:

            | N
        ____|____
        |   v   |
        |  [A]  |
        |   |   |
      N |   v   | N
    ----|->[+]--|--->
    <-- |_______|

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainAdditionNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Gaussian,
                            msg_in1::Message{Gaussian},
                            msg_in2::Any,
                            msg_out::Message{Gaussian})

    dist_1 = ensureParameters!(msg_in1.payload, (:m, :V))
    dist_3 = ensureParameters!(msg_out.payload, (:m, :V))

    outbound_dist.m  = dist_3.m - node.gain[1,1]*dist_1.m
    outbound_dist.V  = dist_3.V + node.gain[1,1]^2 * dist_1.V
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GainAdditionNode:

          ^ | N
        __|_|____
        |   v   |
        |  [A]  |
        |   |   |
      N |   v   | N
    ----|->[+]--|--->
        |_______|

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::GainAdditionNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            msg_in1::Void,
                            msg_in2::Message{Gaussian},
                            msg_out::Message{Gaussian})

    dist_out = ensureParameters!(msg_out.payload, (:m, :V))
    dist_in2 = ensureParameters!(msg_in2.payload, (:m, :V))

    dist_temp = Gaussian(   m = dist_out.m - dist_in2.m,
                                        V = msg_in2.payload.V + msg_out.payload.V)
    ensureParameters!(dist_temp, (:xi, :W))

    outbound_dist.m = NaN
    outbound_dist.V = NaN
    outbound_dist.xi = node.gain[1,1] * dist_temp.xi
    outbound_dist.W = (node.gain[1,1])^2 * dist_temp.W

    return outbound_dist
end


############################################
# MvGaussian methods
############################################

function sumProductRule!{dims_n, dims_m}(   node::GainAdditionNode,
                                            outbound_interface_index::Type{Val{3}},
                                            outbound_dist::MvGaussian{dims_n},
                                            msg_in1::Message{MvGaussian{dims_m}},
                                            msg_in2::Message{MvGaussian{dims_n}},
                                            msg_out::Any)

    return gainAdditionForwardRule!(outbound_dist, msg_in1.payload, msg_in2.payload, node.gain)
end

function sumProductRule!{dims_n, dims_m}(   node::GainAdditionNode,
                                            outbound_interface_index::Type{Val{2}},
                                            outbound_dist::MvGaussian{dims_n},
                                            msg_in1::Message{MvGaussian{dims_m}},
                                            msg_in2::Any,
                                            msg_out::Message{MvGaussian{dims_n}})

    return gainAdditionBackwardIn2Rule!(outbound_dist, msg_in1.payload, msg_out.payload, node.gain)
end

function sumProductRule!{dims_n, dims_m}(   node::GainAdditionNode,
                                            outbound_interface_index::Type{Val{1}},
                                            outbound_dist::MvGaussian{dims_m},
                                            msg_in1::Any,
                                            msg_in2::Message{MvGaussian{dims_n}},
                                            msg_out::Message{MvGaussian{dims_n}})

    dist_temp = vague(MvGaussian{size(node.gain, 1)})
    backwardAdditionRule!(dist_temp, msg_in2.payload, msg_out.payload)
    return gainBackwardRule!(outbound_dist, dist_temp, node.gain, (isdefined(node, :gain_inv)) ? node.gain_inv : nothing)
end

function gainAdditionForwardRule!(dist_result::MvGaussian, dist_1::MvGaussian, dist_2::MvGaussian, A::Any)
    # Select parameterization
    # Order is from least to most computationally intensive
    if isValid(dist_1.m) && isValid(dist_1.V) && isValid(dist_2.m) && isValid(dist_2.V)
        dist_result.m  = forwardGainAdditionMRule(A, dist_2.m, dist_1.m)
        dist_result.V  = forwardGainAdditionVRule(A, dist_2.V, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.m) && isValid(dist_1.W) && isValid(dist_2.m) && isValid(dist_2.W)
        dist_result.m  = forwardGainAdditionMRule(A, dist_2.m, dist_1.m)
        invalidate!(dist_result.V)
        dist_result.W  = forwardGainAdditionWRule(A, dist_2.W, dist_1.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.xi) && isValid(dist_1.W) && isValid(dist_2.xi) && isValid(dist_2.W)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W  = forwardGainAdditionWRule(A, dist_2.W, dist_1.W)
        dist_result.xi = forwardGainAdditionXiRule(A, dist_2.xi, dist_1.xi, dist_2.W, dist_1.W)
    elseif (isValid(dist_1.m) && isValid(dist_1.V)) || (isValid(dist_2.m) && isValid(dist_2.V))
        # Fallback: at least one inbound msg is in (m,V) parametrization
        # Convert the other one to (m,V)
        ensureParameters!(dist_1, (:m, :V))
        ensureParameters!(dist_2, (:m, :V))
        dist_result.m  = forwardGainAdditionMRule(A, dist_2.m, dist_1.m)
        dist_result.V  = forwardGainAdditionVRule(A, dist_2.V, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif (isValid(dist_1.m) && isValid(dist_1.W)) || (isValid(dist_2.m) && isValid(dist_2.W))
        # Fallback: at least one inbound msg is in (m,W) parametrization
        # Convert the other one to (m,W)
        ensureParameters!(dist_1, (:m, :W))
        ensureParameters!(dist_2, (:m, :W))
        dist_result.m  = forwardGainAdditionMRule(A, dist_2.m, dist_1.m)
        invalidate!(dist_result.V)
        dist_result.W  = forwardGainAdditionWRule(A, dist_2.W, dist_1.W)
        invalidate!(dist_result.xi)
    else
        # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
        ensureParameters!(dist_1, (:m, :V))
        ensureParameters!(dist_2, (:m, :V))
        dist_result.m  = forwardGainAdditionMRule(A, dist_2.m, dist_1.m)
        dist_result.V  = forwardGainAdditionVRule(A, dist_2.V, dist_1.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    end

    return dist_result
end

function gainAdditionBackwardIn2Rule!(dist_result::MvGaussian, dist_1::MvGaussian, dist_3::MvGaussian, A::Any)
    # Select parameterization
    # Order is from least to most computationally intensive
    if isValid(dist_1.m) && isValid(dist_1.V) && isValid(dist_3.m) && isValid(dist_3.V)
        dist_result.m  = backwardIn2GainAdditionMRule(A, dist_1.m, dist_3.m)
        dist_result.V  = backwardIn2GainAdditionVRule(A, dist_1.V, dist_3.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.m) && isValid(dist_1.W) && isValid(dist_3.m) && isValid(dist_3.W)
        dist_result.m  = backwardIn2GainAdditionMRule(A, dist_1.m, dist_3.m)
        invalidate!(dist_result.V)
        dist_result.W  = backwardIn2GainAdditionWRule(A, dist_1.W, dist_3.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.xi) && isValid(dist_1.W) && isValid(dist_3.xi) && isValid(dist_3.W)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W  = backwardIn2GainAdditionWRule(A, dist_1.W, dist_3.W)
        dist_result.xi = backwardIn2GainAdditionXiRule(A, dist_1.xi, dist_3.xi, dist_1.W, dist_3.W)
    elseif (isValid(dist_1.m) && isValid(dist_1.V)) || (isValid(dist_3.m) && isValid(dist_3.V))
        # Fallback: at least one inbound msg is in (m,V) parametrization
        # Convert the other one to (m,V)
        ensureParameters!(dist_1, (:m, :V))
        ensureParameters!(dist_3, (:m, :V))
        dist_result.m  = backwardIn2GainAdditionMRule(A, dist_1.m, dist_3.m)
        dist_result.V  = backwardIn2GainAdditionVRule(A, dist_1.V, dist_3.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif (isValid(dist_1.m) && isValid(dist_1.W)) || (isValid(dist_3.m) && isValid(dist_3.W))
        # Fallback: at least one inbound msg is in (m,W) parametrization
        # Convert the other one to (m,W)
        ensureParameters!(dist_1, (:m, :W))
        ensureParameters!(dist_3, (:m, :W))
        dist_result.m  = backwardIn2GainAdditionMRule(A, dist_1.m, dist_3.m)
        invalidate!(dist_result.V)
        dist_result.W  = backwardIn2GainAdditionWRule(A, dist_1.W, dist_3.W)
        invalidate!(dist_result.xi)
    else
        # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
        ensureParameters!(dist_1, (:m, :V))
        ensureParameters!(dist_3, (:m, :V))
        dist_result.m  = backwardIn2GainAdditionMRule(A, dist_1.m, dist_3.m)
        dist_result.V  = backwardIn2GainAdditionVRule(A, dist_1.V, dist_3.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    end

    return dist_result
end

# Rule set for forward propagation ({in1,in2}-->out)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainAdditionMRule{T<:Number}(A::AbstractMatrix{T}, m_x::Vector{T}, m_y::Vector{T}) = m_x + A*m_y
forwardGainAdditionVRule{T<:Number}(A::AbstractMatrix{T}, V_x::AbstractMatrix{T}, V_y::AbstractMatrix{T}) = V_x + A*V_y*A'
forwardGainAdditionWRule{T<:Number}(A::AbstractMatrix{T}, W_x::AbstractMatrix{T}, W_y::AbstractMatrix{T}) = W_x - W_x * A * cholinv(W_y+A'*W_x*A) * A' * W_x
forwardGainAdditionXiRule{T<:Number}(A::AbstractMatrix{T}, xi_x::Vector{T}, xi_y::Vector{T}, W_x::AbstractMatrix{T}, W_y::AbstractMatrix{T}) = xi_x + W_x*A*cholinv(W_y+A'*W_x*A)*(xi_y-A'*xi_x)

# Rule set for backward propagation ({in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardIn2GainAdditionMRule{T<:Number}(A::AbstractMatrix{T}, m_y::Vector{T}, m_z::Vector{T}) = m_z - A*m_y
backwardIn2GainAdditionVRule{T<:Number}(A::AbstractMatrix{T}, V_y::AbstractMatrix{T}, V_z::AbstractMatrix{T}) = V_z + A*V_y*A'
backwardIn2GainAdditionWRule{T<:Number}(A::AbstractMatrix{T}, W_y::AbstractMatrix{T}, W_z::AbstractMatrix{T}) = W_z - W_z * A * cholinv(W_y+A'*W_z*A) * A' * W_z
backwardIn2GainAdditionXiRule{T<:Number}(A::AbstractMatrix{T}, xi_y::Vector{T}, xi_z::Vector{T}, W_y::AbstractMatrix{T}, W_z::AbstractMatrix{T}) = xi_z - W_z*A*cholinv(W_y+A'*W_z*A)*(xi_y+A'*xi_z)
