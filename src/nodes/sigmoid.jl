############################################
# SigmoidNode
############################################
# Description:
#   Links a continuous, real-valued variable (X) to a binary (boolean) one (Y).
#
#     X       Y
#   -----[σ]-----
#
#   f(X,Y) = σ(X⋅Y)
#
# Interfaces:
#   1 i[:real], 2 i[:bool]
#
# Construction:
#   SigmoidNode(:normal_cdf, name="my_node")
#   The optional first argument specifies the sigmoid function.
#   Currently the only option is :normal_cdf.
#
############################################

export SigmoidNode


####################
# Sigmoid functions
####################

# Cummulative Gaussian (CDF of standard normal distribution)
Φ(x::Union(Float64, Vector{Float64})) = 0.5*erfc(-x./sqrt(2.))


####################
# SigmoidNode
####################

type SigmoidNode <: Node
    sigmoid_func::Symbol
    name::ASCIIString
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function SigmoidNode(sigmoid_func::Symbol = :normal_cdf; name = unnamedStr())
        (sigmoid_func == :normal_cdf) || error(":normal_cdf is the only supported sigmoid function at the moment")
        
        self = new(sigmoid_func, name, Array(Interface, 2), Dict{Symbol,Interface}())
        self.i[:real] = self.interfaces[1] = Interface(self)
        self.i[:bool] = self.interfaces[2] = Interface(self)

        return self
    end
end

####################
# Forward messages
####################

function sumProduct!{T<:Real}(node::SigmoidNode,
                              outbound_interface_id::Int,
                              msg_1::Message{DeltaDistribution{T}},
                              ::Nothing)
    # Generate Bernoulli message from incoming Delta message.
    (outbound_interface_id == 2) || error("Invalid call")
    dist_1 = msg_1.payload
    #(typeof(dist_1.m) <: Real) || error("The inbound message should be defined over the 1 dimensional real domain")
    dist_2 = ensureMessage!(node.interfaces[2], BernoulliDistribution).payload

    if node.sigmoid_func == :normal_cdf
        dist_2.p = Φ(dist_1.m)
    else
        error("Unsupported sigmoid function")
    end

    return (:sigmoid_delta_forward,
            node.interfaces[2].message)
end

function sumProduct!(node::SigmoidNode,
                     outbound_interface_id::Int,
                     msg_1::Message{GaussianDistribution},
                     ::Nothing)
    # Generate Bernoulli message from incoming Gaussian message.
    (outbound_interface_id == 2) || error("Invalid call")
    dist_1 = ensureMVParametrization!(msg_1.payload)
    (length(dist_1.m) == 1) || error("Only univariate messages are supported")
    dist_2 = ensureMessage!(node.interfaces[2], BernoulliDistribution).payload

    if node.sigmoid_func == :normal_cdf
        dist_2.p = Φ(dist_1.m[1]/sqrt(1+dist_1.V[1,1]))
    else
        error("Unsupported sigmoid function")
    end

    return (:sigmoid_gaussian_forward,
            node.interfaces[2].message)
end

############################################################
# Backward messages (marginal approximating sum product)
############################################################

function sumProductApprox!{T<:Bool}(node::SigmoidNode,
                                    outbound_interface_id::Int,
                                    msg_1::Message{GaussianDistribution},
                                    msg_2::Message{DeltaDistribution{T}})
    # Generate marginal approximating Gaussian outbound message on i[:real]
    # Procedure:
    #  1. Calculate exact (non-Gaussian) message towards i[:real]
    #  2. Combine exact outbound msg on i[:real] with exact inbound msg on i[:real] to find exact marginal
    #  3. Approximate the exact marginal with a Gaussian one
    #  4. Calculate back the Gaussian outbound msg on i[:real] that yields this approximate Gaussian marginal. Return this msg.
    # NOTE: this calculation results in an implicit cycle in the factor graph since the approximate outbound sumproduct message depends on the inbound message.

    (outbound_interface_id == 1) || error("Invalid call")
    (node.sigmoid_func == :normal_cdf) || error("Unsupported sigmoid function") 

    # Parameters of Gaussian forward message (prior) N(x|μ,σ2)
    dist_1 = ensureMVParametrization!(msg_1.payload)
    μ = dist_1.m[1]; σ2 = dist_1.V[1,1]

    # Parameters of exact backward message Φ(v⋅x)
    v = (msg_2.payload.m==false) ? -1 : 1
    
    # Find the parameters of a Gaussian distribution q(x) that approximates the 'true' marginal p(x)
    # We find the Gaussian by moment matching, as derived in Rasmussen et al. section 3.9.
    z = μ / (v*sqrt(1 + σ2))
    φ(z::Real) = exp(-0.5*z^2)./sqrt(2*pi)
    N = φ(z)
    C = Φ(z)
    q_m = μ + σ2*N / (C*v*sqrt(1+σ2)) # Rasmussen eqn. 3.85
    q_v = σ2 - (σ2^2*N/((1+σ2)*C)) * (z + N/C) # Rasmussen eqn. 3.87 
    node.interfaces[1].edge.marginal = GaussianDistribution(m=q_m, V=q_v)

    # Calculate the approximate message towards i[:real]
    dist_backward = ensureMessage!(node.interfaces[1], GaussianDistribution).payload
    invalidate!(dist_backward.m); invalidate!(dist_backward.V)
    ensureXiWParametrization!(node.interfaces[1].edge.marginal)
    ensureXiWParametrization!(dist_1)
    dist_backward.W = node.interfaces[1].edge.marginal.W - dist_1.W
    dist_backward.xi = node.interfaces[1].edge.marginal.xi - dist_1.xi
 
    return (:empty,
            node.interfaces[1].message)
end

# TODO: sumProductApprox! for BernoulliDistribution instead of DeltaDistribution