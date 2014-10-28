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
#   Also the GaussianNode can accept a fixed mean. This removes the mean interface
#   and stores the mean in the node itself.
#
############################################
#
#   GaussianNode with variable mean:
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
#       GaussianNode(name="my_node")
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
#       Float64 (marginal)
#       Message{Float64}
#
#   Sending:
#   1. (mean):
#       Message{GaussianDistribution}
#       Message{StudentsTDistribution}
#   2. (precision / variance):
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   3. (out):
#       Message{GaussianDistribution}
#
############################################
#
# GaussianNode with fixed mean:
#
#             out
#   ----->[N]----->
#  precision/
#  variance
#
#   out = Message(GaussianDistribution(m=mean, W=precision))
#
#   Example:
#       GaussianNode(name="my_node", m=1.0)
#
# Interface ids, (names) and supported message types:
#   Receiving:
#   1. (precision / variance):
#       Message{Float64}
#       GammaDistribution (marginal)
#       InverseGammaDistribution (marginal)
#   2. (out):
#       Message{Float64}
#
#   Sending:
#   1. (precision / variance):
#       Message{GammaDistribution}
#       Message{InverseGammaDistribution}
#   2. (out):
#       Message{GaussianDistribution}
#
############################################
#
# GaussianNode with fixed variance:
#
#    mean     out
#   ----->[N]----->
#
#   out = Message(GaussianDistribution(m=mean, V=variance))
#
#   Example:
#       GaussianNode(name="my_node", variance=1.0)
#
# Interface ids, (names) and supported message types:
#   Receiving:
#   1. (mean):
#       GaussianDistribution
#   2. (out):
#       GaussianDistribution
#
#   Sending:
#   1. (mean):
#       Message{GaussianDistribution}
#   2. (out):
#       Message{GaussianDistribution}
#
############################################

export GaussianNode

type GaussianNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}

    # Helper fields filled by constructor
    mean::Interface
    precision::Interface
    variance::Interface
    out::Interface
    # Fixed parameters
    m::Vector{Float64}
    V::Matrix{Float64}

    function GaussianNode(; name="unnamed", form::ASCIIString="moment", m::Union(Float64,Vector{Float64},Nothing)=nothing, V::Union(Float64,Matrix{Float64},Nothing)=nothing)
        if m!=nothing && V!=nothing
            error("Only one interface (mean or variance) may be fixed.")
        elseif m!=nothing || V!=nothing
            total_interfaces = 2
        else
            total_interfaces = 3
        end
        self = new(name, Array(Interface, total_interfaces))
        next_interface_index = 1 # Counter keeping track of constructed interfaces

        if m != nothing
            # GaussianNode with fixed mean
            self.m = (typeof(m)==Float64) ? [m] : deepcopy(m)
        else
            # GaussianNode with variable mean
            self.interfaces[next_interface_index] = Interface(self) # Mean interface
            self.mean = self.interfaces[next_interface_index]
            next_interface_index += 1
        end

        # Pick a form for the variance/precision
        if V != nothing
            # GaussianNode with fixed variance
            self.V = (typeof(V)==Float64) ? reshape([V], 1, 1) : deepcopy(V)
        else
            # GaussianNode with variable variance
            self.interfaces[next_interface_index] = Interface(self)
            if form == "moment"
                # Parameters m, V
                self.variance = self.interfaces[next_interface_index]
            elseif form == "precision"
                # Parameters m, W
                self.precision = self.interfaces[next_interface_index]
            elseif form == "canonical"
                error("Canonical form not implemented")
            else
                error("Unrecognized form, $(form). Please use \"moment\", \"canonical\", or \"precision\"")
            end
            next_interface_index += 1
        end
        
        # Out interface
        self.interfaces[next_interface_index] = Interface(self)
        self.out = self.interfaces[next_interface_index]

        return self
    end
end


############################################
# Standard update functions
############################################

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            msg_mean::Any,
                            msg_variance::Message{InverseGammaDistribution},
                            msg_out::Any)
    # Forward / point estimate / InverseGamma
    #                                    
    #            IG         <--      IG  
    #  ---->[N]<----       ---->[N]<---- 
    #  Flt   |               N   |       
    #      N | |                 | Flt     
    #        v v                 v      

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Formulas from table 5.2 in Korl (2005)
    if is(node.interfaces[outbound_interface_id], node.out)
        # Forward message
        gamma = msg_variance.payload
        mean = msg_mean.payload
        dist_out.m = [mean]
        dist_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    elseif is(node.interfaces[outbound_interface_id], node.mean)
        # Backward over mean edge
        # Rules not in Korl, but equivalent by symmetry
        gamma = msg_variance.payload
        y = msg_out.payload
        dist_out.m = [y]
        dist_out.V = reshape([gamma.b/(gamma.a-1)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{InverseGammaDistribution},
                            msg_mean::Message{Float64},
                            ::Nothing,
                            msg_out::Message{Float64})
    # Backward over variance / point estimate / InverseGamma
    #
    #   Flt      IG            
    #  ---->[N]<----
    #        |   -->
    #     Flt|  
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.variance)
        # Backward over variance edge
        y = msg_out.payload
        m = msg_mean.payload
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            msg_variance::Message{InverseGammaDistribution},
                            msg_out::Nothing)
    # Forward with fixed mean
    #
    # IG     -->
    # --->[N]--->
    #         N
    
    # Rule from Korl table 5.2

    (length(node.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    dist_out.m = deepcopy(node.m)
    dist_out.xi = nothing
    dist_out.V = reshape([msg_variance.payload.b/(msg_variance.payload.a-1)], 1, 1)
    dist_out.W = nothing

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{InverseGammaDistribution},
                            msg_variance::Nothing,
                            msg_out::Message{Float64})
    # Backward with fixed mean
    #
    # <--    Flt
    # --->[N]--->
    #  IG

    # Rule from Korl table 5.2

    (length(node.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    dist_out.a = -0.5
    dist_out.b = 0.5*(msg_out.payload - node.m[1])^2

    return node.interfaces[outbound_interface_id].message
end

############################################
# Naive variational update functions
############################################

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Any,
                            marg_variance::InverseGammaDistribution,
                            marg_out::Float64)
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   <--     Q~IG            
    #  ---->[N]<----
    #    N   |   
    #        |Q~Flt  
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.mean) # Mean estimation from variance and sample
        y_0 = marg_out # observation
        a = marg_variance.a # gamma message
        b = marg_variance.b
        dist_out.m = [y_0]
        dist_out.V = reshape([((a+1)/b)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            marg_out::GaussianDistribution)
    # Backward variational update function with fixed variance
    #
    #   <--     Q~N            
    #  ---->[N]---->
    #    N
    #

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.mean)
        ensureMDefined!(marg_out)
        dist_out.m = deepcopy(marg_out.m)
        dist_out.V = deepcopy(node.V)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            marg_precision::GammaDistribution,
                            marg_out::Float64)
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   <--     Q~Gam            
    #  ---->[N]<----
    #    N   |   
    #        |Q~Flt 
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.mean) # Mean estimation from variance and sample
        y_0 = marg_out # observation
        a = marg_precision.a # gamma distribution
        b = marg_precision.b
        dist_out.m = [y_0]
        dist_out.V = nothing
        dist_out.xi = nothing
        dist_out.W = reshape([a/b], 1, 1)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GammaDistribution},
                            marg_mean::GaussianDistribution,
                            ::Nothing,
                            marg_out::Float64)
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N      Gam            
    #  ---->[N]<----
    #       |   -->
    #  Q~Flt|  
    #       v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.precision) # Precision estimation from mean and sample
        y_0 = marg_out # observation
        ensureMWParametrization!(marg_mean)
        (length(marg_mean.W) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        m = marg_mean.m[1] # Gaussian distribution
        W = marg_mean.W[1,1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y_0-m)^2+0.5*inv(W)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GammaDistribution},
                            marg_mean::GaussianDistribution,
                            ::Nothing,
                            marg_out::GaussianDistribution)
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N      Gam            
    #  ---->[N]<----
    #        |  -->
    #    Q~N |  
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.precision) # Precision estimation from mean and sample
        ensureMVParametrization!(marg_mean)
        ensureMVParametrization!(marg_out)
        (length(marg_mean.V) == 1 && length(marg_out.V) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        m = marg_mean.m[1] # Gaussian distribution
        V_m = marg_mean.V[1,1]
        y = marg_out.m[1] # Gaussian distribution
        V_y = marg_out.V[1,1]
        dist_out.a = 1.5
        dist_out.b = 0.5*((y-m)^2+V_m+V_y)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{InverseGammaDistribution},
                            marg_mean::GaussianDistribution,
                            ::Nothing,
                            marg_out::Float64)
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N      IG            
    #  ---->[N]<----
    #        |   -->
    #   Q~Flt|  
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.variance) # Variance estimation from mean and sample
        y_0 = marg_out # observation
        ensureMVParametrization!(marg_mean)
        (length(marg_mean.V) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        m = marg_mean.m[1] # Gaussian message
        V = marg_mean.V[1,1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y_0-m)^2+0.5*V
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            marg_mean::GaussianDistribution,
                            marg_var::InverseGammaDistribution,
                            ::Nothing)
    # Forward over out, InverseGamma input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N     Q~IG            
    #  ---->[N]<----
    #        |
    #      N | |  
    #        v v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.out)
        ensureMDefined!(marg_mean)
        (length(marg_mean.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_mean.m)
        dist_out.V = reshape([marg_var.b/(marg_var.a-1.0)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            marg_var::InverseGammaDistribution,
                            marg_y::GaussianDistribution)
    # Backward over mean, InverseGamma input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # By symmetry the rules are the same as the forward over out.
    #
    #   N       Q~IG            
    #  ---->[N]<----
    #   <--  |
    #    Q~N |  
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.mean)
        ensureMDefined!(marg_y)
        (length(marg_y.W) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_y.m)
        dist_out.V = reshape([marg_var.b/(marg_var.a-1.0)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            marg_mean::GaussianDistribution,
                            ::Nothing)
    # Backward variational update function with fixed variance
    #
    #   Q~N     -->            
    #  ---->[N]---->
    #            N
    #

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.out)
        ensureMDefined!(marg_mean)
        dist_out.m = deepcopy(marg_mean.m)
        dist_out.V = deepcopy(node.V)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            marg_mean::GaussianDistribution,
                            marg_prec::GammaDistribution,
                            ::Nothing)
    # Forward over out, Gamma input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N     Q~Gam            
    #  ---->[N]<----
    #        |
    #      N | |  
    #        v v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.out)
        ensureMDefined!(marg_mean)
        (length(marg_mean.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_mean.m)
        dist_out.W = reshape([marg_prec.a/marg_prec.b], 1, 1)
        dist_out.xi = nothing
        dist_out.V = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            marg_prec::GammaDistribution,
                            marg_y::GaussianDistribution)
    # Backward over mean, Gamma input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # By symmetry the rules are the same as the forward over out.
    #
    #   N       Q~Gam            
    #  ---->[N]<----
    #   <--  |
    #    Q~N |  
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.mean)
        ensureMDefined!(marg_y)
        (length(marg_y.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_y.m)
        dist_out.W = reshape([marg_prec.a/marg_prec.b], 1, 1)
        dist_out.xi = nothing
        dist_out.V = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end


############################################
# Structured variational update functions
############################################

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{StudentsTDistribution},
                            ::Nothing,
                            msg_prec::Message{GammaDistribution},
                            marg_out::GaussianDistribution)
    # Backward message over the mean
    #
    #   <--     Gam            
    #  ---->[N]<----
    #   St   |
    # - - - - - - - -   
    #        | Q~N
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.mean)
        ensureMWParametrization!(marg_out)
        length(marg_out.m) == 1 || error("SVMP update on GaussianNode only defined for univariate distributions")
        m_y = marg_out.m[1]
        prec_y = marg_out.W[1, 1]
        a = msg_prec.payload.a
        b = msg_prec.payload.b
        dist_out.m = [m_y]
        dist_out.W = reshape([(2.0*prec_y*a)/(2.0*prec_y*b + 1)], 1, 1)
        dist_out.nu = 2.0*a
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GammaDistribution},
                            msg_mean::Message{GaussianDistribution},
                            ::Nothing,
                            marg_out::GaussianDistribution)
    # Backward message over the precision
    #
    #    N      Gam            
    #  ---->[N]<----
    #        |   -->
    # - - - - - - - -   
    #        | Q~N
    #        v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.precision)
        ensureMWParametrization!(marg_out)
        length(marg_out.m) == 1 || error("SVMP update on GaussianNode only defined for univariate distributions")
        m_y = marg_out.m[1]
        prec_y = marg_out.W[1, 1]
        ensureMDefined!(msg_mean.payload)
        m_m = msg_mean.payload.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*((1.0/prec_y) + (m_m - m_y)^2)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::GaussianNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            marg::NormalGammaDistribution,
                            ::NormalGammaDistribution, # Same distribution as marg
                            ::Nothing)
    # Forward message over out.
    # This update function has two argments instead of three because it uses the node's joint marginal over the mean and precision.
    #
    #      Q~NGam            
    #  ---->[N]<----
    #        |
    # - - - - - - - -   
    #      | | N  
    #      v v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    if is(node.interfaces[outbound_interface_id], node.out)
        dist_out.m = [marg.m]
        dist_out.W = reshape([marg.a/marg.b], 1, 1)
        dist_out.xi = nothing
        dist_out.V = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end