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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GaussianDistribution, msg_1::Any, msg_2::Message{GaussianDistribution}, msg_3::Message{GaussianDistribution}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Any, msg_3::Message{GaussianDistribution}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Message{GaussianDistribution}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


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
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{T}, msg_1::Any, msg_2::Message{DeltaDistribution{T}}, msg_3::Message{DeltaDistribution{T}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{T}, msg_1::Message{DeltaDistribution{T}}, msg_2::Any, msg_3::Message{DeltaDistribution{T}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{T}, msg_1::Message{DeltaDistribution{T}}, msg_2::Message{DeltaDistribution{T}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::InverseGammaDistribution, msg_1::Any, msg_2::Message{InverseGammaDistribution}, msg_3::Message{InverseGammaDistribution}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::InverseGammaDistribution, msg_1::Message{InverseGammaDistribution}, msg_2::Any, msg_3::Message{InverseGammaDistribution}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::InverseGammaDistribution, msg_1::Message{InverseGammaDistribution}, msg_2::Message{InverseGammaDistribution}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GammaDistribution, msg_1::Any, msg_2::Message{GammaDistribution}, msg_3::Message{GammaDistribution}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GammaDistribution, msg_1::Message{GammaDistribution}, msg_2::Any, msg_3::Message{GammaDistribution}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GammaDistribution, msg_1::Message{GammaDistribution}, msg_2::Message{GammaDistribution}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::BetaDistribution, msg_1::Any, msg_2::Message{BetaDistribution}, msg_3::Message{BetaDistribution}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::BetaDistribution, msg_1::Message{BetaDistribution}, msg_2::Any, msg_3::Message{BetaDistribution}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::BetaDistribution, msg_1::Message{BetaDistribution}, msg_2::Message{BetaDistribution}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::BernoulliDistribution, msg_1::Any, msg_2::Message{BernoulliDistribution}, msg_3::Message{BernoulliDistribution}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::BernoulliDistribution, msg_1::Message{BernoulliDistribution}, msg_2::Any, msg_3::Message{BernoulliDistribution}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::BernoulliDistribution, msg_1::Message{BernoulliDistribution}, msg_2::Message{BernoulliDistribution}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# MvGaussianDistribution methods
############################################

sumProductRule!{T<:MvGaussianDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvGaussianDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvGaussianDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# MvDeltaDistribution methods
############################################

sumProductRule!{T<:MvDeltaDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvDeltaDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvDeltaDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)

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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GaussianDistribution, msg_1::Any, msg_2::Message{GaussianDistribution}, msg_3::Message{StudentsTDistribution}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::GaussianDistribution, msg_1::Any, msg_2::Message{StudentsTDistribution}, msg_3::Message{GaussianDistribution}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Any, msg_3::Message{StudentsTDistribution}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::GaussianDistribution, msg_1::Message{StudentsTDistribution}, msg_2::Any, msg_3::Message{GaussianDistribution}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GaussianDistribution, msg_1::Message{GaussianDistribution}, msg_2::Message{StudentsTDistribution}, msg_3::Any, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::GaussianDistribution, msg_1::Message{StudentsTDistribution}, msg_2::Message{GaussianDistribution}, msg_3::Any, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_2.payload, msg_1.payload, outbound_dist)

function studentsGaussianProd!(dist_gauss_in::GaussianDistribution, dist_stud_in::StudentsTDistribution, dist_result::GaussianDistribution)
    # The student's t distribution is approximated with a Gaussian by moment matching.
    # The result is a Gaussian approximation to the exact result.

    ensureParameters!(dist_gauss_in, (:xi, :W))
    if 0.0 < dist_stud_in.nu <= 1.0
        # The mean and variance for the Student's t are undefined for nu <= 1.
        # However, since we apply a gaussian approximation we assume variance is huge in this case,
        # so the influence of the student's t distribution is negligible.
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

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{GaussianDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{GaussianDistribution}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GaussianDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{GaussianDistribution}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GaussianDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{GaussianDistribution}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# MvGaussian-MvDeltaDistribution combination
############################################

sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MvDeltaDistribution{Float64, dims}, msg_1::Any, msg_2::Message{MvGaussianDistribution{dims}}, msg_3::Message{MvDeltaDistribution{Float64, dims}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MvDeltaDistribution{Float64, dims}, msg_1::Any, msg_2::Message{MvDeltaDistribution{Float64, dims}}, msg_3::Message{MvGaussianDistribution{dims}}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MvDeltaDistribution{Float64, dims}, msg_1::Message{MvGaussianDistribution{dims}}, msg_2::Any, msg_3::Message{MvDeltaDistribution{Float64, dims}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MvDeltaDistribution{Float64, dims}, msg_1::Message{MvDeltaDistribution{Float64, dims}}, msg_2::Any, msg_3::Message{MvGaussianDistribution{dims}}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MvDeltaDistribution{Float64, dims}, msg_1::Message{MvGaussianDistribution{dims}}, msg_2::Message{MvDeltaDistribution{Float64, dims}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MvDeltaDistribution{Float64, dims}, msg_1::Message{MvDeltaDistribution{Float64, dims}}, msg_2::Message{MvGaussianDistribution{dims}}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


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
sumProductRule!{T<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:WishartDistribution}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


#############################################
# Wishart-MatrixDeltaDistribution combination
#############################################

sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MatrixDeltaDistribution{Float64, dims, dims}, msg_1::Any, msg_2::Message{WishartDistribution{dims}}, msg_3::Message{MatrixDeltaDistribution{Float64, dims, dims}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MatrixDeltaDistribution{Float64, dims, dims}, msg_1::Any, msg_2::Message{MatrixDeltaDistribution{Float64, dims, dims}}, msg_3::Message{WishartDistribution{dims}}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MatrixDeltaDistribution{Float64, dims, dims}, msg_1::Message{WishartDistribution{dims}}, msg_2::Any, msg_3::Message{MatrixDeltaDistribution{Float64, dims, dims}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MatrixDeltaDistribution{Float64, dims, dims}, msg_1::Message{MatrixDeltaDistribution{Float64, dims, dims}}, msg_2::Any, msg_3::Message{WishartDistribution{dims}}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MatrixDeltaDistribution{Float64, dims, dims}, msg_1::Message{WishartDistribution{dims}}, msg_2::Message{MatrixDeltaDistribution{Float64, dims, dims}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MatrixDeltaDistribution{Float64, dims, dims}, msg_1::Message{MatrixDeltaDistribution{Float64, dims, dims}}, msg_2::Message{WishartDistribution{dims}}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# LogNormal-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{LogNormalDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{LogNormalDistribution}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{LogNormalDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{LogNormalDistribution}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{LogNormalDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{LogNormalDistribution}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# Gamma-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{GammaDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{GammaDistribution}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GammaDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{GammaDistribution}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{GammaDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{GammaDistribution}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# InverseGamma-DeltaDistribution combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{InverseGammaDistribution}, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Any, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Message{InverseGammaDistribution}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{InverseGammaDistribution}, msg_2::Any, msg_3::Message{DeltaDistribution{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Any, msg_3::Message{InverseGammaDistribution}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{InverseGammaDistribution}, msg_2::Message{DeltaDistribution{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::DeltaDistribution{Float64}, msg_1::Message{DeltaDistribution{Float64}}, msg_2::Message{InverseGammaDistribution}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)
