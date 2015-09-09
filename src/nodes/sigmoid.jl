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
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function SigmoidNode(sigmoid_func::Symbol=:normal_cdf; id=generateNodeId(SigmoidNode))
        (sigmoid_func == :normal_cdf) || error(":normal_cdf is the only supported sigmoid function at the moment")
        
        self = new(sigmoid_func, id, Array(Interface, 2), Dict{Symbol,Interface}())
        self.i[:real] = self.interfaces[1] = Interface(self)
        self.i[:bool] = self.interfaces[2] = Interface(self)

        addNode!(current_graph, self)
        
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

    return (:sigmoid_delta_forward, node.interfaces[2].message)
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

    return (:sigmoid_gaussian_forward, node.interfaces[2].message)
end

############################################################
# Backward messages (expectation)
############################################################

function backwardGaussianExpectationRule!(  node::SigmoidNode,
                                            msg_context::Message{GaussianDistribution}, 
                                            p::Float64)
    # Write back and return approximate backward message (Gaussian) on i[:real] for Bernoulli[p] input.
    # The approximate bw message encodes an 'expectation' under the context of the incoming message N(.|μ,σ2).
    # Procedure:
    #  1. Calculate exact (non-Gaussian) message towards i[:real]
    #  2. Combine exact outbound msg on i[:real] with exact inbound msg ('context') to find exact marginal
    #  3. Approximate the exact (non Gaussian) marginal with a Gaussian one
    #  4. Calculate back the Gaussian outbound msg on i[:real] that yields this approximate Gaussian marginal.
    # NOTE: this calculation results in an implicit cycle in the factor graph since the outbound expectation message depends on the inbound message (context).    
    
    (node.sigmoid_func == :normal_cdf) || error("Unsupported sigmoid function") 
    (0<=p<=1) || error("Bernoulli parameter p should be ∈ [0,1].")

    # Collect parameters of incoming messages (Gaussian and Bernoulli)
    dist_context = ensureMVParametrization!(msg_context.payload)
    μ = dist_context.m[1]; σ2 = dist_context.V[1,1]
    
    # Calculate first and second moment (m_1, m_2) of the 'true' marginal p(x) on edge connected to i[:real]
    # p(x) = f(x) / m_0
    # f(x) = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(x)*N(x|μ,σ2)
    #      = (1-p)*N(x|μ,σ2) + (2p-1)*Φ(z)*(Φ(x)*N(x|μ,σ2)/Φ(z))
    #      = (1-p)*f1(x)     + (2p-1)*Φ(z)*f2(x)
    # See paper for detailed derivation

    # Shorthand symbols
    z = μ / sqrt(1 + σ2)
    N = exp(-0.5*z^2)./sqrt(2*pi) # N(z|0,1)
    C = Φ(z)

    # Moments of f2 = 1/C * Φ(x)*N(x|μ,σ2)
    m2_1 = μ + σ2*N / (C*sqrt(1+σ2))  # First moment of f2 (Rasmussen eqn. 3.85)
    m2_2 = 2*μ*m2_1 - μ^2 + σ2 - (σ2^2*z*N) / (C*sqrt(1+σ2))  # Second moment of f2 (Rasmussen eqn. 3.86)

    # Moments of f(x) (exact marginal)
    m_0 = 1 - p + (2*p-1)*C
    m_1 = ((1-p)*μ + (2*p-1)*C*m2_1)/m_0
    m_2 = ((1-p)*(μ^2+σ2) + (2*p-1)*C*m2_2)/m_0

    # Save Gaussian marginal with the same moments (moment matching)
    marginal = ensureMarginal!(node.interfaces[1].edge, GaussianDistribution)
    invalidate!(marginal.W); invalidate!(marginal.xi)
    marginal.m[1] = m_1
    marginal.V[1,1] = m_2 - m_1^2

    # Calculate the approximate message towards i[:real]
    dist_backward = ensureMessage!(node.interfaces[1], GaussianDistribution).payload
    invalidate!(dist_backward.m); invalidate!(dist_backward.V)
    ensureXiWParametrization!(marginal)
    ensureXiWParametrization!(dist_context)
    dist_backward.W = marginal.W - dist_context.W
    dist_backward.xi = marginal.xi - dist_context.xi
 
    return (:sigmoid_backward_gaussian_expectation, dist_backward)
end

function expectation!{T<:Bool}(node::SigmoidNode,
                                outbound_interface_id::Int,
                                msg_1::Message{GaussianDistribution},
                                msg_2::Message{DeltaDistribution{T}})
    (outbound_interface_id == 1) || error("Invalid call")

    p = (msg_2.payload.m==false) ? 0. : 1.
 
    return backwardGaussianExpectationRule!(node, msg_1, p)
end

function expectation!(node::SigmoidNode,
                      outbound_interface_id::Int,
                      msg_1::Message{GaussianDistribution},
                      msg_2::Message{BernoulliDistribution})
    (outbound_interface_id == 1) || error("Invalid call")

    return backwardGaussianExpectationRule!(node, msg_1, msg_2.payload.p)
end