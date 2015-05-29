############################################
# GainEqualityNode
############################################
# Description:
#   Gain-equality node: A⁻¹*out = in1 = in2
#   Combines the node functions of the FixedGainNode
#   and the EqualityNode for computational efficiency.
#
#        _________
#    in1 |       | in2
#   -----|->[=]<-|-----
#        |   |   |
#        |   v   |
#        |  [A]  |
#        |___|___|
#            | out
#            v
#
#   f(in1,in2,out) = δ(A*in1 - out)⋅δ(A*in2 - out)
#
# Construction:
#   GainEqualityNode([1.0], name="my_node")
#
############################################

export GainEqualityNode

type GainEqualityNode <: Node
    A::Array{Float64}
    name::ASCIIString
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible

    function GainEqualityNode(A::Union(Array{Float64},Float64)=1.0; name=unnamedStr())
        self = new(ensureMatrix(deepcopy(A)), name, Array(Interface, 3), Dict{Symbol,Interface}())

        for (iface_id, iface_name) in enumerate([:in1, :in2, :out])
            self.i[iface_name] = self.interfaces[iface_id] = Interface(self)
        end

        try
            self.A_inv = inv(self.A)
        catch
            warn("The specified multiplier for $(typeof(self)) $(self.name) is not invertible. This might cause problems. Double check that this is what you want.")
        end

        return self
    end
end

isDeterministic(::GainEqualityNode) = true

############################################
# GaussianDistribution methods
############################################

# Rule set for backward propagation ({in2,out}-->in1 or {in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y
backwardGainEqualityVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x - V_x * A' * inv(V_y + A * V_x * A') * A * V_x
backwardGainEqualityMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, V_x::Array{T, 2}, m_y::Array{T, 1}, V_y::Array{T, 2}) = m_x + V_x * A' * inv(V_y + A * V_x * A') * (m_y - A * m_x)

function sumProduct!(node::GainEqualityNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Nothing)
    # Forward message (towards out)
    (outbound_interface_id == 3) || error("The outbound interface id does not match with the calling signature.")
    dist_temp = GaussianDistribution()
    equalityGaussianRule!(dist_temp, msg_in1.payload, msg_in2.payload)
    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload
    fixedGainGaussianForwardRule!(dist_out, dist_temp, node.A, (isdefined(node, :A_inv)) ? node.A_inv : nothing)

    return (:gain_equality_gaussian_forward,
            node.interfaces[outbound_interface_id].message)
end

function sumProduct!(node::GainEqualityNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Nothing,
                            msg_out::Message{GaussianDistribution})
    # Backward message (towards in2)
    return applyBackwardRule!(node, outbound_interface_id, msg_in1, msg_out)
end

function sumProduct!(node::GainEqualityNode,
                            outbound_interface_id::Int,
                            msg_in1::Nothing,
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Message{GaussianDistribution})
    # Backward message (towards in1)
    return applyBackwardRule!(node, outbound_interface_id, msg_in2, msg_out)
end

function applyBackwardRule!(node::GainEqualityNode,
                            outbound_interface_id::Int,
                            msg_in::Message{GaussianDistribution},
                            msg_out::Message{GaussianDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    # Backward message (towards in1 or in2)
    dist_result = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload
    dist_3 = msg_out.payload
    dist_in = msg_in.payload

    # Select parameterization
    # Order is from least to most computationally intensive
    if isValid(dist_3.xi) && isValid(dist_3.W) && isValid(dist_in.xi) && isValid(dist_in.W)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
        dist_result.xi = backwardGainEqualityXiRule(node.A, dist_in.xi, dist_3.xi)
    elseif isValid(dist_3.m) && isValid(dist_3.V) && isValid(dist_in.m) && isValid(dist_in.V)
        dist_result.m = backwardGainEqualityMRule(node.A, dist_in.m, dist_in.V, dist_3.m, dist_3.V)
        dist_result.V = backwardGainEqualityVRule(node.A, dist_in.V, dist_3.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_3.m) && isValid(dist_3.W) && isValid(dist_in.m) && isValid(dist_in.W)
        dist_result.m = backwardGainEqualityMRule(node.A, dist_in.m, inv(dist_in.W), dist_3.m, inv(dist_3.W))
        invalidate!(dist_result.V)
        dist_result.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
        invalidate!(dist_result.xi)
    else
        # Fallback: convert inbound messages to (xi,W) parametrization and then use efficient rules
        ensureXiWParametrization!(dist_in)
        ensureXiWParametrization!(dist_3)
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
        dist_result.xi = backwardGainEqualityXiRule(node.A, dist_in.xi, dist_3.xi)
    end

    return (:gain_equality_gaussian_backward,
            node.interfaces[outbound_interface_id].message)
end