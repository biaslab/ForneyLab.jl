export EqualityNode

"""
Description:

    Equality constraint on 3 variables: i[1] = i[2] = i[3]

         i[2]
         |
    i[1]  |  i[3]
    ------[=]-----

    f(i1,i2,i3) = δ(i1-i3)⋅δ(i2-i3)

Interfaces:

    1. i[1], 2. i[2], 3. i[3]

Construction:

    EqualityNode(id=:my_node)
"""
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

"""
EqualityNode:

      N     N 
    --->[=]<---
         | | 
       N v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GaussianDistribution, msg_1::Any, msg_2::Message{GaussianDistribution}, msg_3::Message{GaussianDistribution}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Any, msg_3::Message{GaussianDistribution}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Message{GaussianDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

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


############################################
# DeltaDistribution methods
############################################

"""
EqualityNode:

      δ     δ 
    --->[=]<---
         | | 
       δ v v
"""
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{T}, msg_1::Any, msg_2::Message{DeltaDistribution{T}}, msg_3::Message{DeltaDistribution{T}}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{T}, msg_1::Message{DeltaDistribution{T}}, msg_2::Any, msg_3::Message{DeltaDistribution{T}}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{T}, msg_1::Message{DeltaDistribution{T}}, msg_2::Message{DeltaDistribution{T}}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

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


############################################
# InverseGammaDistribution methods
############################################

"""
EqualityNode:

     Ig     Ig 
    --->[=]<---
         | | 
      Ig v v

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::InverseGammaDistribution, msg_1::Any, msg_2::Message{InverseGammaDistribution}, msg_3::Message{InverseGammaDistribution}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::InverseGammaDistribution, msg_1::Message{InverseGammaDistribution}, msg_2::Any, msg_3::Message{InverseGammaDistribution}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::InverseGammaDistribution, msg_1::Message{InverseGammaDistribution}, msg_2::Message{InverseGammaDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

function equalityRule!(dist_result::InverseGammaDistribution, dist_1::InverseGammaDistribution, dist_2::InverseGammaDistribution)
    dist_result.a = dist_1.a+dist_2.a+1.0
    dist_result.b = dist_1.b+dist_2.b

    return dist_result
end


############################################
# GammaDistribution methods
############################################

"""
EqualityNode:

    Gam     Gam 
    --->[=]<---
         | | 
     Gam v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GammaDistribution, msg_1::Any, msg_2::Message{GammaDistribution}, msg_3::Message{GammaDistribution}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GammaDistribution, msg_1::Message{GammaDistribution}, msg_2::Any, msg_3::Message{GammaDistribution}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GammaDistribution, msg_1::Message{GammaDistribution}, msg_2::Message{GammaDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

function equalityRule!(dist_result::GammaDistribution, dist_1::GammaDistribution, dist_2::GammaDistribution)
    dist_result.a = dist_1.a+dist_2.a-1.0
    dist_result.b = dist_1.b+dist_2.b

    return dist_result
end


############################################
# BetaDistribution methods
############################################

"""
EqualityNode:

    Beta     Beta 
     --->[=]<---
          | | 
     Beta v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::BetaDistribution, msg_1::Any, msg_2::Message{BetaDistribution}, msg_3::Message{BetaDistribution}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::BetaDistribution, msg_1::Message{BetaDistribution}, msg_2::Any, msg_3::Message{BetaDistribution}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::BetaDistribution, msg_1::Message{BetaDistribution}, msg_2::Message{BetaDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

function equalityRule!(dist_result::BetaDistribution, dist_1::BetaDistribution, dist_2::BetaDistribution)
    dist_result.a = dist_1.a+dist_2.a-1.0
    dist_result.b = dist_1.b+dist_2.b-1.0

    return dist_result
end


############################################
# BernoulliDistribution methods
############################################

"""
EqualityNode:

    Bern     Bern 
     --->[=]<---
          | | 
     Bern v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::BernoulliDistribution, msg_1::Any, msg_2::Message{BernoulliDistribution}, msg_3::Message{BernoulliDistribution}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::BernoulliDistribution, msg_1::Message{BernoulliDistribution}, msg_2::Any, msg_3::Message{BernoulliDistribution}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::BernoulliDistribution, msg_1::Message{BernoulliDistribution}, msg_2::Message{BernoulliDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

function equalityRule!(dist_result::BernoulliDistribution, dist_1::BernoulliDistribution, dist_2::BernoulliDistribution)
    # The result of the Bernoulli equality rule applied to dist_1 and dist_2 is written to dist_result
    norm = dist_1.p * dist_2.p + (1 - dist_1.p) * (1 - dist_2.p)
    (norm > 0) || error("equalityRule! for BernoulliDistribution is not wel l-defined (invalid normalization constant)")
    dist_result.p = (dist_1.p * dist_2.p) / norm
    return dist_result
end


############################################
# MvGaussianDistribution methods
############################################

sumProductRule!{T<:MvGaussianDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!{T<:MvGaussianDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!{T<:MvGaussianDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

function equalityRule!(dist_result::MvGaussianDistribution, dist_1::MvGaussianDistribution, dist_2::MvGaussianDistribution)
    # The result of the Gaussian equality rule applied to dist_1 and dist_2 is written to dist_result
    # The following update rules correspond to node 1 from Table 4.1 in:
    # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.

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

# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
equalityMRule{T<:Number}(m_x::Array{T, 1}, m_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = pinv(W_x+W_y)*(W_x*m_x+W_y*m_y)
equalityVRule{T<:Number}(V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x*pinv(V_x+V_y)*V_y
equalityWRule{T<:Number}(W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x+W_y
equalityXiRule{T<:Number}(xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x+xi_y


############################################
# MvDeltaDistribution methods
############################################

sumProductRule!{T<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!{T<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!{T<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

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


############################################
# Gaussian-Student's t combination
############################################

"""
EqualityNode:

      N     St 
    --->[=]<---
         | | 
       N v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GaussianDistribution, msg_1::Any, msg_2::Message{GaussianDistribution}, msg_3::Message{StudentsTDistribution}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GaussianDistribution, msg_1::Any, msg_2::Message{StudentsTDistribution}, msg_3::Message{GaussianDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Any, msg_3::Message{StudentsTDistribution}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GaussianDistribution, msg_1::Message{StudentsTDistribution}, msg_2::Any, msg_3::Message{GaussianDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Message{StudentsTDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GaussianDistribution, msg_1::Message{StudentsTDistribution}, msg_2::Message{GaussianDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::GaussianDistribution, dist_gauss_in::GaussianDistribution, dist_stud_in::StudentsTDistribution)
    # The result of the Gaussian-Student's t equality rule applied to dist_gauss_in and dist_stud_in is written to dist_result
    # The student's t distribution is approximated by a Gaussian by moment matching.
    # The result is a Gaussian approximation to the exact result.

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


############################################
# Gaussian-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{GaussianDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{GaussianDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GaussianDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{GaussianDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GaussianDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{GaussianDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::DeltaDistribution{Float64}, ::GaussianDistribution, dist_delta::DeltaDistribution{Float64})
    dist_result.m = dist_delta.m
    return dist_result
end

function equalityRule!(dist_result::GaussianDistribution, ::GaussianDistribution, dist_delta::DeltaDistribution{Float64})
    dist_result.m = dist_delta.m
    dist_result.V = tiny
    dist_result.xi = NaN
    dist_result.W = NaN

    return dist_result
end


############################################
# MvGaussian-MvDeltaDistribution combination
############################################

sumProductRule!{TG<:MvGaussianDistribution, TD<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::TD, msg_1::Any, msg_2::Message{TG}, msg_3::Message{TD}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!{TG<:MvGaussianDistribution, TD<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::TD, msg_1::Any, msg_2::Message{TD}, msg_3::Message{TG}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!{TG<:MvGaussianDistribution, TD<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::TD, msg_1::Message{TG}, msg_2::Any, msg_3::Message{TD}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!{TG<:MvGaussianDistribution, TD<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::TD, msg_1::Message{TD}, msg_2::Any, msg_3::Message{TG}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!{TG<:MvGaussianDistribution, TD<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::TD, msg_1::Message{TG}, msg_2::Message{TD}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!{TG<:MvGaussianDistribution, TD<:MvDeltaDistribution{Float64}}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::TD, msg_1::Message{TD}, msg_2::Message{TG}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::MvDeltaDistribution{Float64}, ::MvGaussianDistribution, dist_delta::MvDeltaDistribution{Float64})
    dist_result.m = deepcopy(dist_delta.m)
    dist_result.V = tiny*eye(length(dist_delta.m))
    invalidate!(dist_result.xi)
    invalidate!(dist_result.W)

    return dist_result
end

function equalityRule!(dist_result::MvDeltaDistribution{Float64}, ::MvGaussianDistribution, dist_delta::MvDeltaDistribution{Float64})
    dist_result.m = deepcopy(dist_delta.m)
    return dist_result
end


############################################
# WishartDistribution methods
############################################

"""
EqualityNode:

      W     W 
    --->[=]<---
         | | 
       W v v
"""
sumProductRule!{T<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!{T<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!{T<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)

function equalityRule!(dist_result::WishartDistribution, dist_1::WishartDistribution, dist_2::WishartDistribution)
    # Derivation available in notebook
    dist_result.V = inv(inv(dist_1.V) + inv(dist_2.V))
    dist_result.nu = dist_1.nu + dist_2.nu - size(dist_1.V, 1) - 1.0

    return dist_result
end


############################################
# LogNormal-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{LogNormalDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{LogNormalDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{LogNormalDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{LogNormalDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{LogNormalDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{LogNormalDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::DeltaDistribution{Float64}, dist_1::LogNormalDistribution, dist_2::DeltaDistribution{Float64})
    (dist_2.m >= 0) || error("Can not perform equality rule for log-normal-delta combination for negative numbers")
    dist_result.m = dist_2.m
    return dist_result
end


############################################
# Gamma-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{GammaDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{GammaDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GammaDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{GammaDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GammaDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{GammaDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::DeltaDistribution{Float64}, dist_1::GammaDistribution, dist_2::DeltaDistribution{Float64})
    (dist_2.m >= 0) || error("Can not perform equality rule for gamma-delta combination for negative numbers")
    dist_result.m = dist_2.m
    return dist_result
end


############################################
# InverseGamma-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{InverseGammaDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{InverseGammaDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{InverseGammaDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{InverseGammaDistribution}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{InverseGammaDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{InverseGammaDistribution}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::DeltaDistribution{Float64}, dist_1::InverseGammaDistribution, dist_2::DeltaDistribution{Float64})
    (dist_2.m >= 0) || error("Can not perform equality rule for inverse gamma-delta combination for negative numbers")
    dist_result.m = dist_2.m
    return dist_result
end


############################################
# Wishart-MvDeltaDistribution combination
############################################

sumProductRule!{TD<:MvDeltaDistribution{Float64}, TW<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::TD, msg_1::Any, msg_2::Message{TW}, msg_3::Message{TD}) = return equalityRule!(outbound_dist, msg_2.payload, msg_3.payload)
sumProductRule!{TD<:MvDeltaDistribution{Float64}, TW<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::TD, msg_1::Any, msg_2::Message{TD}, msg_3::Message{TW}) = return equalityRule!(outbound_dist, msg_3.payload, msg_2.payload)
sumProductRule!{TD<:MvDeltaDistribution{Float64}, TW<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::TD, msg_1::Message{TW}, msg_2::Any, msg_3::Message{TD}) = return equalityRule!(outbound_dist, msg_1.payload, msg_3.payload)
sumProductRule!{TD<:MvDeltaDistribution{Float64}, TW<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::TD, msg_1::Message{TD}, msg_2::Any, msg_3::Message{TW}) = return equalityRule!(outbound_dist, msg_3.payload, msg_1.payload)
sumProductRule!{TD<:MvDeltaDistribution{Float64}, TW<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::TD, msg_1::Message{TW}, msg_2::Message{TD}, msg_3::Any) = return equalityRule!(outbound_dist, msg_1.payload, msg_2.payload)
sumProductRule!{TD<:MvDeltaDistribution{Float64}, TW<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::TD, msg_1::Message{TD}, msg_2::Message{TW}, msg_3::Any) = return equalityRule!(outbound_dist, msg_2.payload, msg_1.payload)

function equalityRule!(dist_result::MvDeltaDistribution{Float64}, dist_1::WishartDistribution, dist_2::MvDeltaDistribution{Float64})
    (sum(dist_2.m .< 0) == 0) || error("Can not perform equality rule for wishart-delta combination for negative numbers")
    (size(dist_1.V, 1) == length(dist_2.m)) || error("Dimensions for Wishart and delta must agree")
    dist_result.m = dist_2.m
    return dist_result
end
