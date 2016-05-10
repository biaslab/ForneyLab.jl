export AdditionNode

"""
Description:

    out = in1 + in2

           in2
           |
     in1   v  out
    ----->[+]----->

    f(in1,in2,out) = δ(out - in1 - in2)

Interfaces:

    1 i[:in1], 2 i[:in2], 3 i[:out]

Construction:

    AdditionNode(id=:my_node)
"""
type AdditionNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function AdditionNode(; id=generateNodeId(AdditionNode))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:in1, :in2, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::AdditionNode) = true


############################################
# Gaussian methods
############################################

"""
AdditionNode:

     x~N   y~N 
    --->[+]<---
         | | z~N  
         v v    

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::AdditionNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            msg_in1::Message{Gaussian},
                            msg_in2::Message{Gaussian},
                            msg_out::Any)

    dist_in1 = ensureParameters!(msg_in1.payload, (:m, :V))
    dist_in2 = ensureParameters!(msg_in2.payload, (:m, :V))
    outbound_dist.m = dist_in1.m + dist_in2.m
    outbound_dist.V = dist_in1.V + dist_in2.V
    outbound_dist.W = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end

"""
AdditionNode:

     x~N   y~N 
    --->[+]<---
         |  -->   
     z~N v    

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::AdditionNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Gaussian,
                            msg_in1::Message{Gaussian},
                            msg_in2::Any,
                            msg_out::Message{Gaussian})

    backwardAdditionRule!(outbound_dist, msg_in1.payload, msg_out.payload)
    return outbound_dist
end

"""
AdditionNode:

     x~N   y~N 
    --->[+]<---
    <--  |   
         v z~N    

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 4.1
"""
function sumProductRule!(   node::AdditionNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            msg_in1::Any,
                            msg_in2::Message{Gaussian},
                            msg_out::Message{Gaussian})

    backwardAdditionRule!(outbound_dist, msg_in2.payload, msg_out.payload)
    return outbound_dist
end

function backwardAdditionRule!(dist_result::Gaussian, dist_in::Gaussian, dist_out::Gaussian)
    ensureParameters!(dist_out, (:m, :V))
    ensureParameters!(dist_in, (:m, :V))
    dist_result.m = dist_out.m - dist_in.m
    dist_result.V = dist_out.V + dist_in.V
    dist_result.W = NaN
    dist_result.xi = NaN

    return dist_result
end


#############################################
# Delta methods
#############################################

"""
AdditionNode:

     x~δ   y~δ 
    --->[+]<---
         | |  
     z~δ v v
"""
function sumProductRule!(node::AdditionNode,
                         outbound_interface_index::Type{Val{3}},
                         outbound_dist::Delta{Float64},
                         msg_in1::Message{Delta{Float64}},
                         msg_in2::Message{Delta{Float64}},
                         msg_out::Any)

    outbound_dist.m = msg_in1.payload.m + msg_in2.payload.m
    return outbound_dist
end

"""
AdditionNode:

     x~δ   y~δ 
    --->[+]<---
         |  --> 
     z~δ v   
"""
function sumProductRule!(node::AdditionNode,
                         outbound_interface_index::Type{Val{2}},
                         outbound_dist::Delta{Float64},
                         msg_in1::Message{Delta{Float64}},
                         msg_in2::Any,
                         msg_out::Message{Delta{Float64}})

    backwardAdditionRule!(outbound_dist, msg_in1.payload, msg_out.payload)
    return outbound_dist
end

"""
AdditionNode:

     x~δ   y~δ 
    --->[+]<---
    <--  |  
         v z~δ   
"""
function sumProductRule!(node::AdditionNode,
                         outbound_interface_index::Type{Val{1}},
                         outbound_dist::Delta{Float64},
                         msg_in1::Any,
                         msg_in2::Message{Delta{Float64}},
                         msg_out::Message{Delta{Float64}})

    backwardAdditionRule!(outbound_dist, msg_in2.payload, msg_out.payload)
    return outbound_dist
end

function backwardAdditionRule!(dist_result::Delta{Float64}, dist_in::Delta{Float64}, dist_out::Delta{Float64})
    dist_result.m = dist_out.m - dist_in.m
    return dist_result
end


############################################
# Gaussian-Delta combination
############################################

sumProductRule!(node::AdditionNode,
                outbound_interface_index::Type{Val{3}},
                outbound_dist::Gaussian,
                msg_in1::Message{Delta{Float64}},
                msg_in2::Message{Gaussian},
                msg_out::Any) = sumProductRule!(node, outbound_interface_index, outbound_dist, convert(Message{Gaussian}, msg_in1), msg_in2, msg_out)

sumProductRule!(node::AdditionNode,
                outbound_interface_index::Type{Val{3}},
                outbound_dist::Gaussian,
                msg_in1::Message{Gaussian},
                msg_in2::Message{Delta{Float64}},
                msg_out::Any) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, convert(Message{Gaussian}, msg_in2), msg_out)

sumProductRule!(node::AdditionNode,
                outbound_interface_index::Type{Val{1}},
                outbound_dist::Gaussian,
                msg_in1::Any,
                msg_in2::Message{Delta{Float64}},
                msg_out::Message{Gaussian}) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, convert(Message{Gaussian}, msg_in2), msg_out)

sumProductRule!(node::AdditionNode,
                outbound_interface_index::Type{Val{1}},
                outbound_dist::Gaussian,
                msg_in1::Any,
                msg_in2::Message{Gaussian},
                msg_out::Message{Delta{Float64}}) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, msg_in2, convert(Message{Gaussian}, msg_out))

sumProductRule!(node::AdditionNode,
                outbound_interface_index::Type{Val{2}},
                outbound_dist::Gaussian,
                msg_in1::Message{Delta{Float64}},
                msg_in2::Any,
                msg_out::Message{Gaussian}) = sumProductRule!(node, outbound_interface_index, outbound_dist, convert(Message{Gaussian}, msg_in1), msg_in2, msg_out)

sumProductRule!(node::AdditionNode,
                outbound_interface_index::Type{Val{2}},
                outbound_dist::Gaussian,
                msg_in1::Message{Gaussian},
                msg_in2::Any,
                msg_out::Message{Delta{Float64}}) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, msg_in2, convert(Message{Gaussian}, msg_out))


############################################
# MvGaussian methods
############################################

function sumProductRule!{T<:MvGaussian}(node::AdditionNode,
                                                    outbound_interface_index::Type{Val{3}},
                                                    outbound_dist::T,
                                                    msg_in1::Message{T},
                                                    msg_in2::Message{T},
                                                    msg_out::Any)

    forwardAdditionRule!(outbound_dist, msg_in1.payload, msg_in2.payload)
    return outbound_dist
end

function sumProductRule!{T<:MvGaussian}(node::AdditionNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::T,
                                                    msg_in1::Message{T},
                                                    msg_in2::Any,
                                                    msg_out::Message{T})

    backwardAdditionRule!(outbound_dist, msg_in1.payload, msg_out.payload)
    return outbound_dist
end

function sumProductRule!{T<:MvGaussian}(node::AdditionNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::T,
                                                    msg_in1::Any,
                                                    msg_in2::Message{T},
                                                    msg_out::Message{T})

    backwardAdditionRule!(outbound_dist, msg_in2.payload, msg_out.payload)
    return outbound_dist
end

function forwardAdditionRule!(dist_result::MvGaussian, dist_1::MvGaussian, dist_2::MvGaussian)
    # Calculations for a gaussian message type; Korl (2005), table 4.1
    if isValid(dist_1.m) && isValid(dist_1.V) && isValid(dist_2.m) && isValid(dist_2.V)
        dist_result.m = forwardAdditionMRule(dist_1.m, dist_2.m)
        dist_result.V = forwardAdditionVRule(dist_1.V, dist_2.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.m) && isValid(dist_1.W) && isValid(dist_2.m) && isValid(dist_2.W)
        dist_result.m = forwardAdditionMRule(dist_1.m, dist_2.m)
        invalidate!(dist_result.V)
        dist_result.W = forwardAdditionWRule(dist_1.W, dist_2.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_1.xi) && isValid(dist_1.V) && isValid(dist_2.xi) && isValid(dist_2.V)
        invalidate!(dist_result.m)
        dist_result.V = forwardAdditionVRule(dist_1.V, dist_2.V)
        invalidate!(dist_result.W)
        dist_result.xi= forwardAdditionXiRule(dist_1.V, dist_1.xi, dist_2.V, dist_2.xi)
    else
        # Last resort: calculate (m,V) parametrization for both inbound messages
        ensureParameters!(dist_1, (:m, :V))
        ensureParameters!(dist_2, (:m, :V))
        dist_result.m = forwardAdditionMRule(dist_1.m, dist_2.m)
        dist_result.V = forwardAdditionVRule(dist_1.V, dist_2.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    end

    return dist_result
end

function backwardAdditionRule!(dist_result::MvGaussian, dist_in::MvGaussian, dist_out::MvGaussian)
    # Calculations for a gaussian message type; Korl (2005), table 4.1
    if isValid(dist_in.m) && isValid(dist_in.V) && isValid(dist_out.m) && isValid(dist_out.V)
        dist_result.m = backwardAdditionMRule(dist_in.m, dist_out.m)
        dist_result.V = backwardAdditionVRule(dist_in.V, dist_out.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_in.m) && isValid(dist_in.W) && isValid(dist_out.m) && isValid(dist_out.W)
        dist_result.m = backwardAdditionMRule(dist_in.m, dist_out.m)
        invalidate!(dist_result.V)
        dist_result.W = backwardAdditionWRule(dist_in.W, dist_out.W)
        invalidate!(dist_result.xi)
    elseif isValid(dist_in.xi) && isValid(dist_in.V) && isValid(dist_out.xi) && isValid(dist_out.V)
        invalidate!(dist_result.m)
        dist_result.V = backwardAdditionVRule(dist_in.V, dist_out.V)
        invalidate!(dist_result.W)
        dist_result.xi = backwardAdditionXiRule(dist_in.V, dist_in.xi, dist_out.V, dist_out.xi)
    else
        # Last resort: calculate (m,V) parametrization for both inbound messages
        ensureParameters!(dist_in, (:m, :V))
        ensureParameters!(dist_out, (:m, :V))
        dist_result.m = backwardAdditionMRule(dist_in.m, dist_out.m)
        dist_result.V = backwardAdditionVRule(dist_in.V, dist_out.V)
        invalidate!(dist_result.W)
        invalidate!(dist_result.xi)
    end

    return dist_result
end

# Rule set for forward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardAdditionMRule{T<:Number}(m_x::Vector{T}, m_y::Vector{T}) = m_x + m_y
forwardAdditionVRule{T<:Number}(V_x::AbstractMatrix{T}, V_y::AbstractMatrix{T}) = V_x + V_y
forwardAdditionWRule{T<:Number}(W_x::AbstractMatrix{T}, W_y::AbstractMatrix{T}) = W_x * cholinv(W_x + W_y) * W_y
forwardAdditionXiRule{T<:Number}(V_x::AbstractMatrix{T}, xi_x::Vector{T}, V_y::AbstractMatrix{T}, xi_y::Vector{T}) = cholinv(V_x + V_y) * (V_x*xi_x + V_y*xi_y)

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
# The backward propagation merely negates the mean of the present input message (edge X) and uses the same rules to determine the missing input (edge Y)
# For the sake of clarity there is some redundancy between forward and backward rules.
backwardAdditionMRule{T<:Number}(m_x::Vector{T}, m_z::Vector{T}) = m_z - m_x
backwardAdditionVRule{T<:Number}(V_x::AbstractMatrix{T}, V_z::AbstractMatrix{T}) = V_x + V_z
backwardAdditionWRule{T<:Number}(W_x::AbstractMatrix{T}, W_z::AbstractMatrix{T}) = W_x * cholinv(W_x + W_z) * W_z
backwardAdditionXiRule{T<:Number}(V_x::AbstractMatrix{T}, xi_x::Vector{T}, V_z::AbstractMatrix{T}, xi_z::Vector{T}) = cholinv(V_x + V_z) * (V_z*xi_z - V_x*xi_x)


#############################################
# MvDelta methods
#############################################

function sumProductRule!{T<:MvDelta{Float64}}(  node::AdditionNode,
                                                            outbound_interface_index::Type{Val{3}},
                                                            outbound_dist::T,
                                                            msg_in1::Message{T},
                                                            msg_in2::Message{T},
                                                            msg_out::Any)

    outbound_dist.m = msg_in1.payload.m + msg_in2.payload.m
    return outbound_dist
end

function sumProductRule!{T<:MvDelta{Float64}}(  node::AdditionNode,
                                                            outbound_interface_index::Type{Val{2}},
                                                            outbound_dist::T,
                                                            msg_in1::Message{T},
                                                            msg_in2::Any,
                                                            msg_out::Message{T})

    backwardAdditionRule!(outbound_dist, msg_in1.payload, msg_out.payload)
    return outbound_dist
end

function sumProductRule!{T<:MvDelta{Float64}}(  node::AdditionNode,
                                                            outbound_interface_index::Type{Val{1}},
                                                            outbound_dist::T,
                                                            msg_in1::Any,
                                                            msg_in2::Message{T},
                                                            msg_out::Message{T})

    backwardAdditionRule!(outbound_dist, msg_in2.payload, msg_out.payload)
    return outbound_dist
end

function backwardAdditionRule!(dist_result::MvDelta{Float64}, dist_in::MvDelta{Float64}, dist_out::MvDelta{Float64})
    dist_result.m = dist_out.m - dist_in.m
    return dist_result
end


############################################
# Gaussian-MvDelta combination
############################################

sumProductRule!{TD<:MvDelta{Float64}, TG<:MvGaussian}(
    node::AdditionNode,
    outbound_interface_index::Type{Val{3}},
    outbound_dist::TG,
    msg_in1::Message{TD},
    msg_in2::Message{TG},
    msg_out::Any) = sumProductRule!(node, outbound_interface_index, outbound_dist, convert(Message{MvGaussian}, msg_in1), msg_in2, msg_out)

sumProductRule!{TD<:MvDelta{Float64}, TG<:MvGaussian}(
    node::AdditionNode,
    outbound_interface_index::Type{Val{3}},
    outbound_dist::TG,
    msg_in1::Message{TG},
    msg_in2::Message{TD},
    msg_out::Any) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, convert(Message{MvGaussian}, msg_in2), msg_out)

sumProductRule!{TD<:MvDelta{Float64}, TG<:MvGaussian}(
    node::AdditionNode,
    outbound_interface_index::Type{Val{1}},
    outbound_dist::TG,
    msg_in1::Any,
    msg_in2::Message{TD},
    msg_out::Message{TG}) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, convert(Message{MvGaussian}, msg_in2), msg_out)

sumProductRule!{TD<:MvDelta{Float64}, TG<:MvGaussian}(
    node::AdditionNode,
    outbound_interface_index::Type{Val{1}},
    outbound_dist::TG,
    msg_in1::Any,
    msg_in2::Message{TG},
    msg_out::Message{TD}) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, msg_in2, convert(Message{MvGaussian}, msg_out))

sumProductRule!{TD<:MvDelta{Float64}, TG<:MvGaussian}(
    node::AdditionNode,
    outbound_interface_index::Type{Val{2}},
    outbound_dist::TG,
    msg_in1::Message{TD},
    msg_in2::Any,
    msg_out::Message{TG}) = sumProductRule!(node, outbound_interface_index, outbound_dist, convert(Message{MvGaussian}, msg_in1), msg_in2, msg_out)

sumProductRule!{TD<:MvDelta{Float64}, TG<:MvGaussian}(
    node::AdditionNode,
    outbound_interface_index::Type{Val{2}},
    outbound_dist::TG,
    msg_in1::Message{TG},
    msg_in2::Any,
    msg_out::Message{TD}) = sumProductRule!(node, outbound_interface_index, outbound_dist, msg_in1, msg_in2, convert(Message{MvGaussian}, msg_out))
