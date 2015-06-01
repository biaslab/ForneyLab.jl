############################################
# GainAdditionNode
############################################
# Description:
#   Gain-addition node: out = A*in1 + in2
#   Combines the node functions of the FixedGainNode
#   and the AdditionNode for computational efficiency.
#
#            | in1
#            |
#        ____|____
#        |   v   |
#        |  [A]  |
#        |   |   |
#    in2 |   v   | out
#   -----|->[+]--|---->
#        |_______|
#
#   f(in1,in2,out) = δ(out - A*in1 - in2)
#
# Interfaces:
#   1 i[:in1], 2 i[:in2], 3 i[:out]
#   
# Construction:
#   GainAdditionNode([1.0], id=:my_node)
#
############################################
export GainAdditionNode

type GainAdditionNode <: Node
    A::Array{Float64}
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible

    function GainAdditionNode(A::Union(Array{Float64},Float64)=1.0; id=generateNodeId())
        self = new(ensureMatrix(deepcopy(A)), id, Array(Interface, 3), Dict{Symbol,Interface}())
        !haskey(current_graph.n, id) ? current_graph.n[id] = self : error("Node id $(id) already present")

        for (iface_id, iface_handle) in enumerate([:in1, :in2, :out])
            self.i[iface_handle] = self.interfaces[iface_id] = Interface(self)
        end

        try
            self.A_inv = inv(self.A)
        catch
            warn("The specified multiplier for $(typeof(self)) $(self.id) is not invertible. This might cause problems. Double check that this is what you want.")
        end

        return self
    end
end

isDeterministic(::GainAdditionNode) = true

############################################
# GaussianDistribution methods
############################################

# Rule set for forward propagation ({in1,in2}-->out)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainAdditionMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, m_y::Array{T, 1}) = m_x + A*m_y
forwardGainAdditionVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x + A*V_y*A'
forwardGainAdditionWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x - W_x * A * inv(W_y+A'*W_x*A) * A' * W_x
forwardGainAdditionXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = xi_x + W_x*A*inv(W_y+A'*W_x*A)*(xi_y-A'*xi_x)

# Rule set for backward propagation ({in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardIn2GainAdditionMRule{T<:Number}(A::Array{T, 2}, m_y::Array{T, 1}, m_z::Array{T, 1}) = m_z - A*m_y
backwardIn2GainAdditionVRule{T<:Number}(A::Array{T, 2}, V_y::Array{T, 2}, V_z::Array{T, 2}) = V_z + A*V_y*A'
backwardIn2GainAdditionWRule{T<:Number}(A::Array{T, 2}, W_y::Array{T, 2}, W_z::Array{T, 2}) = W_z - W_z * A * inv(W_y+A'*W_z*A) * A' * W_z
backwardIn2GainAdditionXiRule{T<:Number}(A::Array{T, 2}, xi_y::Array{T, 1}, xi_z::Array{T, 1}, W_y::Array{T, 2}, W_z::Array{T, 2}) = xi_z - W_z*A*inv(W_y+A'*W_z*A)*(xi_y+A'*xi_z)

# Forward to OUT
function sumProduct!(node::GainAdditionNode,
                            outbound_interface_id::Int,
                            in1::Message{GaussianDistribution},
                            in2::Message{GaussianDistribution},
                            ::Nothing)

    if outbound_interface_id == 3
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_1 = in1.payload
        dist_2 = in2.payload

        # Select parameterization
        # Order is from least to most computationally intensive
        if isValid(dist_1.m) && isValid(dist_1.V) && isValid(dist_2.m) && isValid(dist_2.V)
            dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
            dist_out.V  = forwardGainAdditionVRule(node.A, dist_2.V, dist_1.V)
            invalidate!(dist_out.W) 
            invalidate!(dist_out.xi)
        elseif isValid(dist_1.m) && isValid(dist_1.W) && isValid(dist_2.m) && isValid(dist_2.W)
            dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
            invalidate!(dist_out.V) 
            dist_out.W  = forwardGainAdditionWRule(node.A, dist_2.W, dist_1.W)
            invalidate!(dist_out.xi)
        elseif isValid(dist_1.xi) && isValid(dist_1.W) && isValid(dist_2.xi) && isValid(dist_2.W)
            invalidate!(dist_out.m) 
            invalidate!(dist_out.V) 
            dist_out.W  = forwardGainAdditionWRule(node.A, dist_2.W, dist_1.W)
            dist_out.xi = forwardGainAdditionXiRule(node.A, dist_2.xi, dist_1.xi, dist_2.W, dist_1.W)
        elseif (isValid(dist_1.m) && isValid(dist_1.V)) || (isValid(dist_2.m) && isValid(dist_2.V))
            # Fallback: at least one inbound msg is in (m,V) parametrization
            # Convert the other one to (m,V)
            ensureMVParametrization!(dist_1)
            ensureMVParametrization!(dist_2)
            dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
            dist_out.V  = forwardGainAdditionVRule(node.A, dist_2.V, dist_1.V)
            invalidate!(dist_out.W) 
            invalidate!(dist_out.xi)
        elseif (isValid(dist_1.m) && isValid(dist_1.W)) || (isValid(dist_2.m) && isValid(dist_2.W))
            # Fallback: at least one inbound msg is in (m,W) parametrization
            # Convert the other one to (m,W)
            ensureMWParametrization!(dist_1)
            ensureMWParametrization!(dist_2)
            dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
            invalidate!(dist_out.V) 
            dist_out.W  = forwardGainAdditionWRule(node.A, dist_2.W, dist_1.W)
            invalidate!(dist_out.xi)
        else
            # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
            ensureMVParametrization!(dist_1)
            ensureMVParametrization!(dist_2)
            dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
            dist_out.V  = forwardGainAdditionVRule(node.A, dist_2.V, dist_1.V)
            invalidate!(dist_out.W) 
            invalidate!(dist_out.xi)
        end

        return (:gain_addition_gaussian_forward,
                node.interfaces[outbound_interface_id].message)
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.id).")
    end
end

# Backward to IN2
function sumProduct!(node::GainAdditionNode,
                            outbound_interface_id::Int,
                            in1::Message{GaussianDistribution},
                            ::Nothing,
                            out::Message{GaussianDistribution})

    if outbound_interface_id == 2
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_1 = in1.payload
        dist_3 = out.payload

        # Select parameterization
        # Order is from least to most computationally intensive
        if isValid(dist_1.m) && isValid(dist_1.V) && isValid(dist_3.m) && isValid(dist_3.V)
            dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
            dist_out.V  = backwardIn2GainAdditionVRule(node.A, dist_1.V, dist_3.V)
            invalidate!(dist_out.W) 
            invalidate!(dist_out.xi)
        elseif isValid(dist_1.m) && isValid(dist_1.W) && isValid(dist_3.m) && isValid(dist_3.W)
            dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
            invalidate!(dist_out.V) 
            dist_out.W  = backwardIn2GainAdditionWRule(node.A, dist_1.W, dist_3.W)
            invalidate!(dist_out.xi)
        elseif isValid(dist_1.xi) && isValid(dist_1.W) && isValid(dist_3.xi) && isValid(dist_3.W)
            invalidate!(dist_out.m) 
            invalidate!(dist_out.V) 
            dist_out.W  = backwardIn2GainAdditionWRule(node.A, dist_1.W, dist_3.W)
            dist_out.xi = backwardIn2GainAdditionXiRule(node.A, dist_1.xi, dist_3.xi, dist_1.W, dist_3.W)
        elseif (isValid(dist_1.m) && isValid(dist_1.V)) || (isValid(dist_3.m) && isValid(dist_3.V))
            # Fallback: at least one inbound msg is in (m,V) parametrization
            # Convert the other one to (m,V)
            ensureMVParametrization!(dist_1)
            ensureMVParametrization!(dist_3)
            dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
            dist_out.V  = backwardIn2GainAdditionVRule(node.A, dist_1.V, dist_3.V)
            invalidate!(dist_out.W) 
            invalidate!(dist_out.xi)
        elseif (isValid(dist_1.m) && isValid(dist_1.W)) || (isValid(dist_3.m) && isValid(dist_3.W))
            # Fallback: at least one inbound msg is in (m,W) parametrization
            # Convert the other one to (m,W)
            ensureMWParametrization!(dist_1)
            ensureMWParametrization!(dist_3)
            dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
            invalidate!(dist_out.V) 
            dist_out.W  = backwardIn2GainAdditionWRule(node.A, dist_1.W, dist_3.W)
            invalidate!(dist_out.xi)
        else
            # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
            ensureMVParametrization!(dist_1)
            ensureMVParametrization!(dist_3)
            dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
            dist_out.V  = backwardIn2GainAdditionVRule(node.A, dist_1.V, dist_3.V)
            invalidate!(dist_out.W) 
            invalidate!(dist_out.xi)
        end

        return (:gain_addition_gaussian_backward_in2,
                node.interfaces[outbound_interface_id].message)
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.id).")
    end
end

# Backward to IN1
function sumProduct!(node::GainAdditionNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            in2::Message{GaussianDistribution},
                            out::Message{GaussianDistribution})

    if outbound_interface_id == 1
        dist_temp = GaussianDistribution()
        additionGaussianBackwardRule!(dist_temp, in2.payload, out.payload)
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload
        fixedGainGaussianBackwardRule!(dist_out, dist_temp, node.A, (isdefined(node, :A_inv)) ? node.A_inv : nothing)

        return (:gain_addition_gaussian_backward_in1,
            node.interfaces[outbound_interface_id].message)
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.id).")
    end
end