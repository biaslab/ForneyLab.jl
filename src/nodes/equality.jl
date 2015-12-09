############################################
# EqualityNode
############################################
# Description:
#   Equality constraint on 3 variables: i[1] = i[2] = i[3]
#
#          i[2]
#          |
#    i[1]  |  i[3]
#   ------[=]-----
#
#   f(i1,i2,i3) = δ(i1-i3)⋅δ(i2-i3)
#
# Interfaces:
#   1 i[1], 2 i[2], 3 i[3]
#
# Construction:
#   EqualityNode(id=:my_node)
#
############################################

export EqualityNode

type EqualityNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Int,Interface}

    function EqualityNode(; id=generateNodeId(EqualityNode))
        self = new(id, Array(Interface, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for interface_index = 1:3
            self.i[interface_index] = self.interfaces[interface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::EqualityNode) = true

# Implement firstFreeInterface since EqualityNode is symmetrical in its interfaces
function firstFreeInterface(node::EqualityNode)
    # Return the first free interface of a symmetrical node
    for interface in node.interfaces
        if interface.partner == nothing
            return interface
        end
    end
    error("No free interface on $(typeof(node)) $(node.id)")
end

############################################
# GaussianDistribution methods
############################################

function equalityRule!(dist_result::GaussianDistribution, dist_1::GaussianDistribution, dist_2::GaussianDistribution)
    # The result of the Gaussian equality rule applied to dist_1 and dist_2 is written to dist_result
    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
    ensureParameters!(dist_1, (:xi, :W))
    ensureParameters!(dist_2, (:xi, :W))
    dist_result.m = NaN
    dist_result.V = NaN
    dist_result.W  = dist_1.W + dist_2.W
    dist_result.xi = dist_1.xi + dist_2.xi

    return dist_result
end

function sumProduct!(   node::EqualityNode,
                        outbound_interface_index::Int,
                        msg_1::Message{GaussianDistribution},
                        msg_2::Message{GaussianDistribution},
                        ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.
    dist_result = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    equalityRule!(dist_result, msg_1.payload, msg_2.payload)

    return (:equality_gaussian,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{GaussianDistribution},
            msg_2::Message{GaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{GaussianDistribution},
            ::Void,
            msg_2::Message{GaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

############################################
# Gaussian-Student's t combination
############################################

function equalityRule!(dist_result::GaussianDistribution, dist_gauss_in::GaussianDistribution, dist_stud_in::StudentsTDistribution)
    # The result of the Gaussian-Student's t equality rule applied to dist_gauss_in and dist_stud_in is written to dist_result
    # The student's t distribution is approximated by a Gaussian by moment matching.
    # The result is a Gaussian approximation to the exact result.
    (isProper(dist_gauss_in) && isProper(dist_stud_in)) || error("Inputs of equalityRule! should be proper distributions")

    ensureParameters!(dist_gauss_in, (:xi, :W))
    if 0.0 < dist_stud_in.nu <= 1.0
        # The mean and variance for the Student's t are undefined for nu <= 1.
        # However, since we apply a gaussian approximation we assume variance is huge in this case,
        # so the second incoming message dominates.
        approx_V = huge
        approx_m = dist_stud_in.m
    else
        approx_V = var(dist_stud_in)
        approx_m = mean(dist_stud_in)
    end

    approx_W = inv(approx_V)
    approx_xi = approx_W * approx_m
    dist_result.xi  = dist_gauss_in.xi + approx_xi
    dist_result.W  = dist_gauss_in.W + approx_W
    dist_result.V = NaN
    dist_result.m = NaN

    return dist_result
end

equalityRule!(  dist_result::GaussianDistribution,
                dist_stud_in::StudentsTDistribution,
                dist_gauss_in::GaussianDistribution) = equalityRule!(dist_result, dist_guass_in, dist_stud_in)

function sumProduct!(   node::EqualityNode,
                        outbound_interface_index::Int,
                        msg_st::Message{StudentsTDistribution},
                        msg_n::Message{GaussianDistribution},
                        ::Void)
    # Combination of Gaussian and student's t-distribution
    # Definitions available in derivations notebook
    # Same as Gaussian equality rule

    dist_result = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    equalityRule!(dist_result, msg_n.payload, msg_st.payload)

    return (:equality_gaussian_student,
            node.interfaces[outbound_interface_index].message)
end

# Call signatures for different message orders
# Outbound interface 3
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_n::Message{GaussianDistribution},
            msg_st::Message{StudentsTDistribution},
            ::Void) = sumProduct!(node, outbound_interface_index, msg_st, msg_n, nothing)
# Outbound interface 2
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_n::Message{GaussianDistribution},
            ::Void,
            msg_st::Message{StudentsTDistribution}) = sumProduct!(node, outbound_interface_index, msg_st, msg_n, nothing)
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_st::Message{StudentsTDistribution},
            ::Void,
            msg_n::Message{GaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_st, msg_n, nothing)
# Outbound interface 1
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_n::Message{GaussianDistribution},
            msg_st::Message{StudentsTDistribution}) = sumProduct!(node, outbound_interface_index, msg_st, msg_n, nothing)
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_st::Message{StudentsTDistribution},
            msg_n::Message{GaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_st, msg_n, nothing)

############################################
# DeltaDistribution methods
############################################

function equalityRule!(dist_result::DeltaDistribution, dist_1::DeltaDistribution, dist_2::DeltaDistribution)
    # The result of the delta equality rule applied to dist_1 and dist_2 is written to dist_result
    # Outbound message is equal to the inbound messages if both inbound messages are equal.
    # Otherwise, the outbound message is zero.
    if dist_1 == dist_2
        dist_result.m = deepcopy(dist_1.m)
    else
        dist_result.m = zero(dist_1.m)
    end

    return dist_result
end

function sumProduct!(node::EqualityNode,
                     outbound_interface_index::Int,
                     msg_1::Message{DeltaDistribution{Float64}},
                     msg_2::Message{DeltaDistribution{Float64}},
                     ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.

    dist_result = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload

    equalityRule!(dist_result, msg_1.payload, msg_2.payload)

    return (:equality_delta,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{DeltaDistribution{Float64}},
            ::Void,
            msg_2::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{DeltaDistribution{Float64}},
            msg_2::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)


############################################
# InverseGammaDistribution methods
############################################

function equalityRule!(dist_result::InverseGammaDistribution, dist_1::InverseGammaDistribution, dist_2::InverseGammaDistribution)
    # The result of the inverse gamma equality rule applied to dist_1 and dist_2 is written to dist_result
    # Definition from Korl table 5.2
    (isProper(dist_1) && isProper(dist_2)) || error("Inputs of equalityRule! should be proper distributions")
    dist_result.a = dist_1.a+dist_2.a+1.0
    dist_result.b = dist_1.b+dist_2.b
    return dist_result
end

function sumProduct!(node::EqualityNode,
                            outbound_interface_index::Int,
                            msg_1::Message{InverseGammaDistribution},
                            msg_2::Message{InverseGammaDistribution},
                            ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], InverseGammaDistribution).payload

    equalityRule!(dist_out, msg_1.payload, msg_2.payload)

    return (:equality_inverse_gamma,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{InverseGammaDistribution},
            ::Void,
            msg_2::Message{InverseGammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{InverseGammaDistribution},
            msg_2::Message{InverseGammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

############################################
# GammaDistribution methods
############################################

function equalityRule!(dist_result::GammaDistribution, dist_1::GammaDistribution, dist_2::GammaDistribution)
    # The result of the gamma equality rule applied to dist_1 and dist_2 is written to dist_result
    # Derivation available in notebook
    (isProper(dist_1) && isProper(dist_2)) || error("Inputs of equalityRule! should be proper distributions")
    dist_result.a = dist_1.a+dist_2.a-1.0
    dist_result.b = dist_1.b+dist_2.b
    return dist_result
end

function sumProduct!(node::EqualityNode,
                            outbound_interface_index::Int,
                            msg_1::Message{GammaDistribution},
                            msg_2::Message{GammaDistribution},
                            ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GammaDistribution).payload

    equalityRule!(dist_out, msg_1.payload, msg_2.payload)

    return (:equality_gamma,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{GammaDistribution},
            ::Void,
            msg_2::Message{GammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{GammaDistribution},
            msg_2::Message{GammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

############################################
# BetaDistribution methods
############################################

function equalityRule!(dist_result::BetaDistribution, dist_1::BetaDistribution, dist_2::BetaDistribution)
    # The result of the beta equality rule applied to dist_1 and dist_2 is written to dist_result
    # Derivation available in notebook
    (isProper(dist_1) && isProper(dist_2)) || error("Inputs of equalityRule! should be proper distributions")
    dist_result.a = dist_1.a+dist_2.a-1.0
    dist_result.b = dist_1.b+dist_2.b-1.0
    return dist_result
end

function sumProduct!(node::EqualityNode,
                            outbound_interface_index::Int,
                            msg_1::Message{BetaDistribution},
                            msg_2::Message{BetaDistribution},
                            ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], BetaDistribution).payload

    equalityRule!(dist_out, msg_1.payload, msg_2.payload)

    return (:equality_beta,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{BetaDistribution},
            ::Void,
            msg_2::Message{BetaDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{BetaDistribution},
            msg_2::Message{BetaDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

############################################
# Gaussian-DeltaDistribution combination
############################################

function equalityRule!(dist_result::DeltaDistribution{Float64}, ::GaussianDistribution, dist_delta::DeltaDistribution{Float64})
    dist_result.m = dist_delta.m
    return dist_result
end

equalityRule!(  dist_result::DeltaDistribution{Float64},
                dist_delta::DeltaDistribution{Float64},
                dist_gauss::GaussianDistribution) = equalityRule!(dist_result, dist_gauss, dist_delta)

function equalityRule!(dist_result::GaussianDistribution, ::GaussianDistribution, dist_delta::DeltaDistribution{Float64})
    dist_result.m = dist_delta.m
    dist_result.V = tiny
    dist_result.W = NaN
    dist_result.xi = NaN
    return dist_result
end

equalityRule!(  dist_result::GaussianDistribution,
                dist_delta::DeltaDistribution{Float64},
                dist_gauss::GaussianDistribution) = equalityRule!(dist_result, dist_gauss, dist_delta)

function sumProduct!(node::EqualityNode,
                     outbound_interface_index::Int,
                     msg_delta::Message{DeltaDistribution{Float64}},
                     msg_n::Message{GaussianDistribution},
                     ::Void)
    # Combination of Gaussian and delta
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload

    equalityRule!(dist_out, msg_n.payload, msg_delta.payload)

    return (:equality_gaussian_delta,
            node.interfaces[outbound_interface_index].message)
end

# Call signatures for messages other ways around
# Outbound interface 3
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_n::Message{GaussianDistribution},
            msg_delta::Message{DeltaDistribution{Float64}},
            ::Void) = sumProduct!(node, outbound_interface_index, msg_delta, msg_n, nothing)
# Outbound interface 2
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_n::Message{GaussianDistribution},
            ::Void,
            msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_n, nothing)
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_delta::Message{DeltaDistribution{Float64}},
            ::Void,
            msg_n::Message{GaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_n, nothing)
# Outbound interface 1
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_n::Message{GaussianDistribution},
            msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_n, nothing)
sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_delta::Message{DeltaDistribution{Float64}},
            msg_n::Message{GaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_n, nothing)

############################################
# MvGaussianDistribution methods
############################################

# Rule set
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
equalityMRule{T<:Number}(m_x::Array{T, 1}, m_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = pinv(W_x+W_y)*(W_x*m_x+W_y*m_y)
equalityVRule{T<:Number}(V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x*pinv(V_x+V_y)*V_y
equalityWRule{T<:Number}(W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x+W_y
equalityXiRule{T<:Number}(xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x+xi_y

function equalityRule!(dist_result::MvGaussianDistribution, dist_1::MvGaussianDistribution, dist_2::MvGaussianDistribution)
    # The result of the Gaussian equality rule applied to dist_1 and dist_2 is written to dist_result
    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
    (isProper(dist_1) && isProper(dist_2)) || error("Inputs of equalityRule! should be proper distributions")
    if isValid(dist_1.m) && isValid(dist_1.W) && isValid(dist_2.m) && isValid(dist_2.W)
        dist_result.m  = equalityMRule(dist_1.m, dist_2.m, dist_1.W, dist_2.W)
        invalidate!(dist_result.V)
        dist_result.W  = equalityWRule(dist_1.W, dist_2.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.xi) && isValid(dist_1.V) && isValid(dist_2.xi) && isValid(dist_2.V)
        invalidate!(dist_result.m)
        dist_result.V  = equalityVRule(dist_1.V, dist_2.V)
        invalidate!(dist_result.W)
        dist_result.xi = equalityXiRule(dist_1.xi, dist_2.xi)
    else
        # Use (xi,W)
        ensureParameters!(dist_1, (:xi, :W))
        ensureParameters!(dist_2, (:xi, :W))
        invalidate!(dist_result.m)
        invalidate!(dist_result.V)
        dist_result.W  = equalityWRule(dist_1.W, dist_2.W)
        dist_result.xi = equalityXiRule(dist_1.xi, dist_2.xi)
    end

    return dist_result
end

function sumProduct!(   node::EqualityNode,
                        outbound_interface_index::Int,
                        msg_1::Message{MvGaussianDistribution},
                        msg_2::Message{MvGaussianDistribution},
                        ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.
    dist_result = ensureMessage!(node.interfaces[outbound_interface_index], MvGaussianDistribution).payload

    equalityRule!(dist_result, msg_1.payload, msg_2.payload)

    return (:equality_gaussian,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{MvGaussianDistribution},
            msg_2::Message{MvGaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{MvGaussianDistribution},
            ::Void,
            msg_2::Message{MvGaussianDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

############################################
# MvDeltaDistribution methods
############################################

function equalityRule!(dist_result::MvDeltaDistribution, dist_1::MvDeltaDistribution, dist_2::MvDeltaDistribution)
    # The result of the delta equality rule applied to dist_1 and dist_2 is written to dist_result
    # Outbound message is equal to the inbound messages if both inbound messages are equal.
    # Otherwise, the outbound message is zero.
    if dist_1 == dist_2
        dist_result.m = deepcopy(dist_1.m)
    else
        dist_result.m = zero(dist_1.m)
    end

    return dist_result
end

function sumProduct!(node::EqualityNode,
                     outbound_interface_index::Int,
                     msg_1::Message{MvDeltaDistribution{Float64}},
                     msg_2::Message{MvDeltaDistribution{Float64}},
                     ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.

    dist_result = ensureMessage!(node.interfaces[outbound_interface_index], MvDeltaDistribution{Float64}).payload

    equalityRule!(dist_result, msg_1.payload, msg_2.payload)

    return (:equality_delta,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{MvDeltaDistribution{Float64}},
            ::Void,
            msg_2::Message{MvDeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{MvDeltaDistribution{Float64}},
            msg_2::Message{MvDeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)


############################################
# LogNormal-DeltaDistribution combination
############################################

function equalityRule!(dist_result::DeltaDistribution{Float64}, dist_1::LogNormalDistribution, dist_2::DeltaDistribution{Float64})
    (dist_2.m >= 0) || error("Can not perform equality rule for log-normal-delta combination for negative numbers")
    dist_result.m = dist_2.m
    return dist_result
end

equalityRule!(  dist_result::DeltaDistribution{Float64},
                dist_1::DeltaDistribution{Float64},
                dist_2::LogNormalDistribution) = equalityRule!(dist_result, dist_2, dist_1)

function sumProduct!(node::EqualityNode,
                     outbound_interface_index::Int,
                     msg_delta::Message{DeltaDistribution{Float64}},
                     msg_logn::Message{LogNormalDistribution},
                     ::Void)
    # Combination of gamma and float
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload

    equalityRule!(dist_out, msg_delta.payload, msg_logn.payload)

    return (:equality_log_normal_delta,
            node.interfaces[outbound_interface_index].message)
end
# Call signature for messages other ways around
# Outbound interface 3
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_logn::Message{LogNormalDistribution}, msg_delta::Message{DeltaDistribution{Float64}}, ::Void) = sumProduct!(node, outbound_interface_index, msg_delta, msg_logn, nothing)
# Outbound interface 2
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_logn::Message{LogNormalDistribution}, ::Void, msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_logn, nothing)
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_delta::Message{DeltaDistribution{Float64}}, ::Void, msg_logn::Message{LogNormalDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_logn, nothing)
# Outbound interface 1
sumProduct!(node::EqualityNode, outbound_interface_index::Int, ::Void, msg_logn::Message{LogNormalDistribution}, msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_logn, nothing)
sumProduct!(node::EqualityNode, outbound_interface_index::Int, ::Void, msg_delta::Message{DeltaDistribution{Float64}}, msg_logn::Message{LogNormalDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_logn, nothing)


############################################
# Gamma-DeltaDistribution combination
############################################

function equalityRule!(dist_result::DeltaDistribution{Float64}, dist_1::GammaDistribution, dist_2::DeltaDistribution{Float64})
    (dist_2.m >= 0) || error("Can not perform equality rule for gamma-delta combination for negative numbers")
    dist_result.m = dist_2.m
    return dist_result
end

equalityRule!(  dist_result::DeltaDistribution{Float64},
                dist_1::DeltaDistribution{Float64},
                dist_2::GammaDistribution) = equalityRule!(dist_result, dist_2, dist_1)

function sumProduct!(node::EqualityNode,
                     outbound_interface_index::Int,
                     msg_delta::Message{DeltaDistribution{Float64}},
                     msg_gam::Message{GammaDistribution},
                     ::Void)
    # Combination of gamma and float
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload

    equalityRule!(dist_out, msg_delta.payload, msg_gam.payload)

    return (:equality_gamma_delta,
            node.interfaces[outbound_interface_index].message)
end
# Call signature for messages other ways around
# Outbound interface 3
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_gam::Message{GammaDistribution}, msg_delta::Message{DeltaDistribution{Float64}}, ::Void) = sumProduct!(node, outbound_interface_index, msg_delta, msg_gam, nothing)
# Outbound interface 2
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_gam::Message{GammaDistribution}, ::Void, msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_gam, nothing)
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_delta::Message{DeltaDistribution{Float64}}, ::Void, msg_gam::Message{GammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_gam, nothing)
# Outbound interface 1
sumProduct!(node::EqualityNode, outbound_interface_index::Int, ::Void, msg_gam::Message{GammaDistribution}, msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_gam, nothing)
sumProduct!(node::EqualityNode, outbound_interface_index::Int, ::Void, msg_delta::Message{DeltaDistribution{Float64}}, msg_gam::Message{GammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_gam, nothing)

############################################
# InverseGamma-DeltaDistribution combination
############################################

function equalityRule!(dist_result::DeltaDistribution{Float64}, dist_1::InverseGammaDistribution, dist_2::DeltaDistribution{Float64})
    (dist_2.m >= 0) || error("Can not perform equality rule for inverse gamma-delta combination for negative numbers")
    dist_result.m = dist_2.m
    return dist_result
end
equalityRule!(  dist_result::DeltaDistribution{Float64},
                dist_1::DeltaDistribution{Float64},
                dist_2::InverseGammaDistribution) = equalityRule!(dist_result, dist_2, dist_1)

function sumProduct!(node::EqualityNode,
                     outbound_interface_index::Int,
                     msg_delta::Message{DeltaDistribution{Float64}},
                     msg_igam::Message{InverseGammaDistribution},
                     ::Void)
    # Combination of gamma and float
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], DeltaDistribution{Float64}).payload

    equalityRule!(dist_out, msg_delta.payload, msg_igam.payload)

    return (:equality_inverse_gamma_delta,
            node.interfaces[outbound_interface_index].message)
end
# Call signature for messages other ways around
# Outbound interface 3
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_igam::Message{InverseGammaDistribution}, msg_delta::Message{DeltaDistribution{Float64}}, ::Void) = sumProduct!(node, outbound_interface_index, msg_delta, msg_igam, nothing)
# Outbound interface 2
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_igam::Message{InverseGammaDistribution}, ::Void, msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_igam, nothing)
sumProduct!(node::EqualityNode, outbound_interface_index::Int, msg_delta::Message{DeltaDistribution{Float64}}, ::Void, msg_igam::Message{InverseGammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_igam, nothing)
# Outbound interface 1
sumProduct!(node::EqualityNode, outbound_interface_index::Int, ::Void, msg_igam::Message{InverseGammaDistribution}, msg_delta::Message{DeltaDistribution{Float64}}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_igam, nothing)
sumProduct!(node::EqualityNode, outbound_interface_index::Int, ::Void, msg_delta::Message{DeltaDistribution{Float64}}, msg_igam::Message{InverseGammaDistribution}) = sumProduct!(node, outbound_interface_index, msg_delta, msg_igam, nothing)

############################################
# BernoulliDistribution methods
############################################

function equalityRule!(dist_result::BernoulliDistribution, dist_1::BernoulliDistribution, dist_2::BernoulliDistribution)
    # The result of the Bernoulli equality rule applied to dist_1 and dist_2 is written to dist_result
    norm = dist_1.p * dist_2.p + (1 - dist_1.p) * (1 - dist_2.p)
    (norm > 0) || error("equalityRule! for BernoulliDistribution is not wel l-defined (invalid normalization constant)")
    dist_result.p = (dist_1.p * dist_2.p) / norm
    return dist_result
end

function sumProduct!(   node::EqualityNode,
                        outbound_interface_index::Int,
                        msg_1::Message{BernoulliDistribution},
                        msg_2::Message{BernoulliDistribution},
                        ::Void)
    # Calculate an outbound message based on the inbound messages and the node function.
    dist_result = ensureMessage!(node.interfaces[outbound_interface_index], BernoulliDistribution).payload

    equalityRule!(dist_result, msg_1.payload, msg_2.payload)

    return (:equality_bernoulli,
            node.interfaces[outbound_interface_index].message)
end

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            ::Void,
            msg_1::Message{BernoulliDistribution},
            msg_2::Message{BernoulliDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)

sumProduct!(node::EqualityNode,
            outbound_interface_index::Int,
            msg_1::Message{BernoulliDistribution},
            ::Void,
            msg_2::Message{BernoulliDistribution}) = sumProduct!(node, outbound_interface_index, msg_1, msg_2, nothing)
