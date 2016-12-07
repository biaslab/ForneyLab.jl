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
# Gaussian methods
############################################

"""
EqualityNode:

      N     N
    --->[=]<---
         | |
       N v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Gaussian, msg_1::Any, msg_2::Message{Gaussian}, msg_3::Message{Gaussian}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Gaussian, msg_1::Message{Gaussian}, msg_2::Any, msg_3::Message{Gaussian}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Gaussian, msg_1::Message{Gaussian}, msg_2::Message{Gaussian}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# Delta methods
############################################

"""
EqualityNode:

      δ     δ
    --->[=]<---
         | |
       δ v v
"""
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{T}, msg_1::Any, msg_2::Message{Delta{T}}, msg_3::Message{Delta{T}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{T}, msg_1::Message{Delta{T}}, msg_2::Any, msg_3::Message{Delta{T}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:Float64}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{T}, msg_1::Message{Delta{T}}, msg_2::Message{Delta{T}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# InverseGamma methods
############################################

"""
EqualityNode:

     Ig     Ig
    --->[=]<---
         | |
      Ig v v

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::InverseGamma, msg_1::Any, msg_2::Message{InverseGamma}, msg_3::Message{InverseGamma}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::InverseGamma, msg_1::Message{InverseGamma}, msg_2::Any, msg_3::Message{InverseGamma}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::InverseGamma, msg_1::Message{InverseGamma}, msg_2::Message{InverseGamma}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# Gamma methods
############################################

"""
EqualityNode:

    Gam     Gam
    --->[=]<---
         | |
     Gam v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Gamma, msg_1::Any, msg_2::Message{Gamma}, msg_3::Message{Gamma}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Gamma, msg_1::Message{Gamma}, msg_2::Any, msg_3::Message{Gamma}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Gamma, msg_1::Message{Gamma}, msg_2::Message{Gamma}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# Beta methods
############################################

"""
EqualityNode:

    Beta     Beta
     --->[=]<---
          | |
     Beta v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Beta, msg_1::Any, msg_2::Message{Beta}, msg_3::Message{Beta}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Beta, msg_1::Message{Beta}, msg_2::Any, msg_3::Message{Beta}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Beta, msg_1::Message{Beta}, msg_2::Message{Beta}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# Bernoulli methods
############################################

"""
EqualityNode:

    Bern     Bern
     --->[=]<---
          | |
     Bern v v
"""
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Bernoulli, msg_1::Any, msg_2::Message{Bernoulli}, msg_3::Message{Bernoulli}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Bernoulli, msg_1::Message{Bernoulli}, msg_2::Any, msg_3::Message{Bernoulli}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Bernoulli, msg_1::Message{Bernoulli}, msg_2::Message{Bernoulli}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# Categorical{k} methods
############################################

"""
EqualityNode:
    Cat       Cat
     --->[=]<---
          | |
      Cat v v
"""
sumProductRule!{k}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Categorical{k}, msg_1::Any, msg_2::Message{Categorical{k}}, msg_3::Message{Categorical{k}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{k}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Categorical{k}, msg_1::Message{Categorical{k}}, msg_2::Any, msg_3::Message{Categorical{k}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{k}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Categorical{k}, msg_1::Message{Categorical{k}}, msg_2::Message{Categorical{k}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)

############################################
# Dirichlet methods
############################################

"""
EqualityNode:
    Dir       Dir
     --->[=]<---
          | |
      Dir v v
"""
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Dirichlet{dims}, msg_1::Any, msg_2::Message{Dirichlet{dims}}, msg_3::Message{Dirichlet{dims}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Dirichlet{dims}, msg_1::Message{Dirichlet{dims}}, msg_2::Any, msg_3::Message{Dirichlet{dims}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Dirichlet{dims}, msg_1::Message{Dirichlet{dims}}, msg_2::Message{Dirichlet{dims}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# MvGaussian methods
############################################

sumProductRule!{T<:MvGaussian}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvGaussian}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvGaussian}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


############################################
# MvDelta methods
############################################

sumProductRule!{T<:MvDelta}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvDelta}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:MvDelta}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)

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
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Gaussian, msg_1::Any, msg_2::Message{Gaussian}, msg_3::Message{StudentsT}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Gaussian, msg_1::Any, msg_2::Message{StudentsT}, msg_3::Message{Gaussian}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Gaussian, msg_1::Message{Gaussian}, msg_2::Any, msg_3::Message{StudentsT}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Gaussian, msg_1::Message{StudentsT}, msg_2::Any, msg_3::Message{Gaussian}, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Gaussian, msg_1::Message{Gaussian}, msg_2::Message{StudentsT}, msg_3::Any, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Gaussian, msg_1::Message{StudentsT}, msg_2::Message{Gaussian}, msg_3::Any, approx::Type{MomentMatching}) = studentsGaussianProd!(msg_2.payload, msg_1.payload, outbound_dist)

function studentsGaussianProd!(dist_gauss_in::Gaussian, dist_stud_in::StudentsT, dist_result::Gaussian)
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
# Gaussian-Delta combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{Gaussian}, msg_3::Message{Delta{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{Delta{Float64}}, msg_3::Message{Gaussian}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{Gaussian}, msg_2::Any, msg_3::Message{Delta{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Any, msg_3::Message{Gaussian}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{Gaussian}, msg_2::Message{Delta{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Message{Gaussian}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# MvGaussian-MvDelta combination
############################################

sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MvDelta{Float64, dims}, msg_1::Any, msg_2::Message{MvGaussian{dims}}, msg_3::Message{MvDelta{Float64, dims}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MvDelta{Float64, dims}, msg_1::Any, msg_2::Message{MvDelta{Float64, dims}}, msg_3::Message{MvGaussian{dims}}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MvDelta{Float64, dims}, msg_1::Message{MvGaussian{dims}}, msg_2::Any, msg_3::Message{MvDelta{Float64, dims}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MvDelta{Float64, dims}, msg_1::Message{MvDelta{Float64, dims}}, msg_2::Any, msg_3::Message{MvGaussian{dims}}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MvDelta{Float64, dims}, msg_1::Message{MvGaussian{dims}}, msg_2::Message{MvDelta{Float64, dims}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MvDelta{Float64, dims}, msg_1::Message{MvDelta{Float64, dims}}, msg_2::Message{MvGaussian{dims}}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# Wishart methods
############################################

"""
EqualityNode:

      W     W
    --->[=]<---
         | |
       W v v
"""
sumProductRule!{T<:Wishart}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::T, msg_1::Any, msg_2::Message{T}, msg_3::Message{T}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:Wishart}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::T, msg_1::Message{T}, msg_2::Any, msg_3::Message{T}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{T<:Wishart}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::T, msg_1::Message{T}, msg_2::Message{T}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)


#############################################
# Wishart-MatrixDelta combination
#############################################

sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MatrixDelta{Float64, dims, dims}, msg_1::Any, msg_2::Message{Wishart{dims}}, msg_3::Message{MatrixDelta{Float64, dims, dims}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::MatrixDelta{Float64, dims, dims}, msg_1::Any, msg_2::Message{MatrixDelta{Float64, dims, dims}}, msg_3::Message{Wishart{dims}}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MatrixDelta{Float64, dims, dims}, msg_1::Message{Wishart{dims}}, msg_2::Any, msg_3::Message{MatrixDelta{Float64, dims, dims}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::MatrixDelta{Float64, dims, dims}, msg_1::Message{MatrixDelta{Float64, dims, dims}}, msg_2::Any, msg_3::Message{Wishart{dims}}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MatrixDelta{Float64, dims, dims}, msg_1::Message{Wishart{dims}}, msg_2::Message{MatrixDelta{Float64, dims, dims}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!{dims}(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::MatrixDelta{Float64, dims, dims}, msg_1::Message{MatrixDelta{Float64, dims, dims}}, msg_2::Message{Wishart{dims}}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# LogNormal-Delta combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{LogNormal}, msg_3::Message{Delta{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{Delta{Float64}}, msg_3::Message{LogNormal}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{LogNormal}, msg_2::Any, msg_3::Message{Delta{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Any, msg_3::Message{LogNormal}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{LogNormal}, msg_2::Message{Delta{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Message{LogNormal}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# Gamma-Delta combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{Gamma}, msg_3::Message{Delta{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{Delta{Float64}}, msg_3::Message{Gamma}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{Gamma}, msg_2::Any, msg_3::Message{Delta{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Any, msg_3::Message{Gamma}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{Gamma}, msg_2::Message{Delta{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Message{Gamma}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)


############################################
# InverseGamma-Delta combination
############################################

sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{InverseGamma}, msg_3::Message{Delta{Float64}}) = prod!(msg_2.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{1}}, outbound_dist::Delta{Float64}, msg_1::Any, msg_2::Message{Delta{Float64}}, msg_3::Message{InverseGamma}) = prod!(msg_3.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{InverseGamma}, msg_2::Any, msg_3::Message{Delta{Float64}}) = prod!(msg_1.payload, msg_3.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{2}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Any, msg_3::Message{InverseGamma}) = prod!(msg_3.payload, msg_1.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{InverseGamma}, msg_2::Message{Delta{Float64}}, msg_3::Any) = prod!(msg_1.payload, msg_2.payload, outbound_dist)
sumProductRule!(node::EqualityNode, outbound_interface_index::Type{Val{3}}, outbound_dist::Delta{Float64}, msg_1::Message{Delta{Float64}}, msg_2::Message{InverseGamma}, msg_3::Any) = prod!(msg_2.payload, msg_1.payload, outbound_dist)
