############################################
# GaussianNode
############################################
# Description:
#   Node converting an input mean and precision
#   to a univariate Gaussian distribution.
#   
#   The GaussianNode has two different names for its second interface,
#   namely precision and variance. These named handles are simply
#   pointers to one and the same interface. 
#
#         mean
#          |
#          v  out
#   ----->[N]----->
#  precision/
#  variance
#
#   out = Message(GaussianDistribution(m=mean, W=precision))
#
#   Example:
#       GaussianNode([1.0], [0.1]; name="my_node")
#
# Interface ids, (names) and supported message types:
#   Receiving:
#   1. (mean):
#       Message{Float64}
#       GaussianDistribution (marginal)
#   2. (precision / variance):
#       Message{Float64}
#       GammaDistribution (marginal)
#       InverseGammaDistribution (marginal)
#   3. (out):
#       Message{Float64}
#
#   Sending:
#   1. (mean):
#       Message{GaussianDistribution}
#   2. (precision / variance):
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   3. (out):
#       Message{GaussianDistribution}
############################################

export GaussianNode

type GaussianNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    variational::Bool
    # Helper fields filled by constructor
    mean::Interface
    precision::Interface
    variance::Interface
    out::Interface

    # TODO: THE VARIATIONAL FLAG SHOULD ONLY BE USED AS A MEANS TO PROPAGATE INFORMATION TO THE INTERFACES
    function GaussianNode(variational; name="unnamed", form::ASCIIString="moment", args...)
        self = new(name, Array(Interface, 3), variational)

        # Look for fixed parameters
        args = Dict(zip(args...)...) # Cast args to dictionary
        # Set m if required
        self.interfaces[1] = Interface(self)
        self.mean = self.interfaces[1]
        if haskey(args, :m)
            Edge(ForneyLab.ClampNode(Message(args[:m])).out, self.mean, typeof(args[:m]))
        end
        # Pick a form for the variance/precision
        self.interfaces[2] = Interface(self)
        if form == "moment"
            # Parameters m, V
            self.variance = self.interfaces[2]
            if haskey(args, :V)
                Edge(ForneyLab.ClampNode(Message(args[:V])).out, self.variance, typeof(args[:V]))
            end
        elseif form == "precision"
            # Parameters m, W
            self.precision = self.interfaces[2]
            if haskey(args, :W)
                Edge(ForneyLab.ClampNode(Message(args[:W])).out, self.precision, typeof(args[:W]))
            end
        elseif form == "canonical"
            error("Canonical form not implemented")
        else
            error("Unrecognized form, $(form). Please use \"moment\", \"canonical\", or \"precision\"")
        end
        self.interfaces[3] = Interface(self) # Set out interface
        self.out = self.interfaces[3]

        return self
    end
end
GaussianNode(; args...) = GaussianNode(false; args...)


############################################
# Update functions
############################################

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GaussianDistribution},
                            msg1::Any,
                            msg2::Message{InverseGammaDistribution},
                            msg3::Any)
    # Forward / point estimate / InverseGamma

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Formulas from table 5.2 in Korl (2005)
    if is(node.interfaces[outbound_interface_id], node.out)
        # Forward message
        gamma = msg2.value
        mean = msg1.value
        dist_out.m = [mean]
        dist_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    elseif is(node.interfaces[outbound_interface_id], node.mean)
        # Backward over mean edge
        # Rules not in Korl, but equivalent by symmetry
        gamma = msg2.value
        y = msg3.value
        dist_out.m = [y]
        dist_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{InverseGammaDistribution},
                            msg1::Message{Float64},
                            ::Any,
                            msg3::Message{Float64})
    # Backward over variance / point estimate / InverseGamma

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if is(node.interfaces[outbound_interface_id], node.variance)
        # Backward over variance edge
        y = msg3.value
        m = msg1.value
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2
    else
        error("Message-type ($(inbound_message_types)) outbound_interface_id ($(outbound_interface_id)) combination not defined for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GaussianDistribution},
                            ::Any,
                            marg2::InverseGammaDistribution,
                            msg3::Message{Float64})
    # Backward over mean / variational / InverseGamma
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if is(node.interfaces[outbound_interface_id], node.mean) # Mean estimation from variance and sample
        y_0 = msg3.value # observation
        a = marg2.a # gamma message
        b = marg2.b
        dist_out.m = [y_0]
        dist_out.V = reshape([((a+1)/b)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GaussianDistribution},
                            ::Any,
                            marg2::GammaDistribution,
                            msg3::Message{Float64})
    # Backward over mean / variational / Gamma
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if is(node.interfaces[outbound_interface_id], node.mean) # Mean estimation from variance and sample
        y_0 = msg3.value # observation
        a = marg2.a # gamma distribution
        b = marg2.b
        dist_out.m = [y_0]
        dist_out.V = nothing
        dist_out.xi = nothing
        dist_out.W = reshape([a/b], 1, 1)
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GammaDistribution},
                            marg1::GaussianDistribution,
                            ::Any,
                            msg3::Message{Float64})
    # Backward over precision / variational / Gamma
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if is(node.interfaces[outbound_interface_id], node.precision) # Precision estimation from mean and sample
        y_0 = msg3.value # observation
        ensureMWParametrization!(marg1)
        m = marg1.m[1] # Gaussian distribution
        W = marg1.W[1,1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y_0-m)^2+0.5*inv(W)
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{InverseGammaDistribution},
                            marg1::GaussianDistribution,
                            ::Any,
                            msg3::Message{Float64})
    # Backward over variance / variational / InverseGamma
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    if is(node.interfaces[outbound_interface_id], node.variance) # Variance estimation from mean and sample
        y_0 = msg3.value # observation
        ensureMVParametrization!(marg1)
        m = marg1.m[1] # Gaussian message
        V = marg1.V[1,1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y_0-m)^2+0.5*V
    else
        error("Inbound message type $(inbound_message_types) outbound_interface_id $(outbound_interface_id) undefined for type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end