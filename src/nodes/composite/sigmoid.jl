############################################
# SigmoidCompositeNode
############################################
# Description:
#   Transformation through a softened sigmoid function
#
#      ----------------
#      |   |N(0, |a|b |
#  in1 |   vgam) v v  | out
#   ---|->[+]-->[ S ]-|--->
#      |              |
#      ----------------
#
#   out = sigmoid(in1 + N(0, gamma); a, b),
#         with sigmoid(x; a, b) = 1/(1 + a*exp(-b*x))
#
#   Example:
#       SigmoidCompositeNode(a, b; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       Message{GaussianDistribution}
#       GaussianDistribution
#   2. (out):
#       Message{BetaDistribution}
#       BetaDistribution
############################################

export SigmoidCompositeNode

type SigmoidCompositeNode <: Node
    use_composite_update_rules::Bool
    a::Float64
    b::Float64
    gamma::Float64
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    out::Interface

    function SigmoidCompositeNode(use_composite_update_rules::Bool=true; a=1.0, b=1.0, gamma=huge(), name=unnamedStr())
        use_composite_update_rules == true || error("SigmoidCompositeNode $(name) does not support explicit internal message passing")
        gamma > 0.0 || error("Gamma for SigmoidCompositeNode $(name) must be positive")
        self = new(use_composite_update_rules, a, b, gamma, name, Array(Interface, 2))

        # Set up the interfaces
        param_list = [:in1, :out]
        for i = 1:length(param_list)
            self.interfaces[i] = Interface(self) # Construct interface
            setfield!(self, param_list[i], self.interfaces[i]) # Set named interfaces
        end

        return self
    end
end

isDeterministic(::SigmoidCompositeNode) = false

############################################
# Standard update functions
############################################

function sumProduct!(node::SigmoidCompositeNode,
                     outbound_interface_id::Int,
                     msg_in1::Message{GaussianDistribution},
                     ::Nothing)

    ensureMWParametrization!(msg_in1.payload)
    (length(msg_in1.payload.m) == 1) || error("SigmoidCompositeNode only implemented for unvariate distributions")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], BetaDistribution).payload
    # Numeric optimization by minimizing KL divergence between beta and the unnormalized analytic outgoing distribution

    # Compute the precision of the message that enters the sigmoid node
    W_mid = forwardAdditionWRule(reshape([node.gamma],1,1), msg_in1.payload.W)

    # Note, the optimization objective KLBetaq is a closure that requires x and q to be set in this parent scope
    x = linspace(0.01,0.99,100) # Note, the resulting answer depends on the resolution of x
    q = exp(-0.5*W_mid[1,1]*(node.b-msg_in1.payload.m[1]-(1/node.a)*(log((1./x)-1)) ).^2) # Unnormalized PDF of analytic outgoing message 

    # Closure wrapping KL divergence calculation for beta and params (a, b)
    KLBetaq(params::Vector{Float64}) = ForneyLab.KLBetaq(x, q; a=params[1], b=params[2])

    optimum = optimize(KLBetaq, [2.0, 2.0])

    (a_est, b_est) = optimum.minimum
    dist_out.a = a_est
    dist_out.b = b_est

    return (:sigmoid_forward,
            node.interfaces[outbound_interface_id].message)
end

function sumProduct!(node::SigmoidCompositeNode,
                     outbound_interface_id::Int,
                     ::Nothing,
                     msg_out::Message{BetaDistribution})

    dist_in1 = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    # See sigmoid derivations notebook, calculated through Laplace aproximation
    # TODO: this Laplace approximation may be too poor
    dist_out=msg_out.payload
    (dist_out.a > 1.0 && dist_out.b > 1.0 && node.a > 0.0) || error("Backward sumproduct update rule for sigmoid node only defined for Beta a and b > 0 and sigmoid a > 0")

    # Calculate the message backwards out of the sigmoid
    W_mid = reshape([ ((dist_out.a-1)/(dist_out.a+dist_out.b-2))^dist_out.a * ((dist_out.b-1)/(dist_out.a+dist_out.b-2))^dist_out.b * (dist_out.a+dist_out.b-2) * node.b^2 ],1,1)
    # Add the noise of gamma to the outcome
    W_in1 = backwardAdditionWRule(reshape([node.gamma],1,1), W_mid)

    dist_in1.m = [(1.0/node.b) * log(node.a*( (dist_out.a-1.0)/(dist_out.b-1.0) ))]
    dist_in1.V = nothing
    dist_in1.xi = nothing
    dist_in1.W = W_in1 

    return (:sigmoid_backward,
            node.interfaces[outbound_interface_id].message)
end


############################################
# Variational update functions
############################################

function sumProduct!(node::SigmoidCompositeNode,
                     outbound_interface_id::Int,
                     dist_in1::GaussianDistribution,
                     ::Nothing)

    ensureMWParametrization!(dist_in1)
    (length(dist_in1.m) == 1) || error("SigmoidCompositeNode only implemented for unvariate distributions")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], BetaDistribution).payload
    # Numeric optimization by minimizing KL divergence between beta and the unnormalized analytic outgoing distribution

    # Note, the optimization objective KLBetaq is a closure that requires x and q to be set in this parent scope
    x = linspace(0.01,0.99,100) # Note, the resulting answer depends on the resolution of x
    q = ( (1./x)-1 ).^( (node.gamma*(node.b+dist_in1.m[1])/node.a) - (node.gamma*log((1./x)-1))/(2.0*node.a^2) ) # Unnormalized PDF of analytic outgoing message 

    # Closure wrapping KL divergence calculation for beta and params (a, b)
    KLBetaq(params::Vector{Float64}) = ForneyLab.KLBetaq(x, q; a=params[1], b=params[2])

    optimum = optimize(KLBetaq, [2.0, 2.0])

    (a_est, b_est) = optimum.minimum
    dist_out.a = a_est
    dist_out.b = b_est

    return (:sigmoid_forward_variational,
            node.interfaces[outbound_interface_id].message)
end

function sumProduct!(node::SigmoidCompositeNode,
                     outbound_interface_id::Int,
                     ::Nothing,
                     dist_out::BetaDistribution)

    dist_in1 = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    # See sigmoid derivations notebook, calculated through Laplace aproximation
    (dist_out.a > 0.0 && dist_out.b > 0.0 && node.a > 0.0 && node.gamma > 0.0) || error("Backward variational update rule for sigmoid node only defined for Beta a and b > 0 and sigmoid a > 0")

    dist_in1.m = [ (1.0/node.b) * log( node.a*(dist_out.a/dist_out.b) ) ]
    dist_in1.V = nothing
    dist_in1.xi = nothing
    dist_in1.W = reshape([ dist_out.a^2*dist_out.b^2*node.gamma*node.b^2 / (dist_out.a+dist_out.b)^4 ],1,1)

    return (:sigmoid_backward_variational,
            node.interfaces[outbound_interface_id].message)
end
