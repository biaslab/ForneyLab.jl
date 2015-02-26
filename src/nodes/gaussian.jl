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
#       Message{DeltaDistribution}
#       GaussianDistribution (marginal)
#   2. (precision / variance):
#       Message{DeltaDistribution}
#       GammaDistribution (marginal)
#       InverseGammaDistribution (marginal)
#   3. (out):
#       Message{DeltaDistribution}
#       GaussianDistribution (marginal)
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
#       Message{DeltaDistribution}
#       GammaDistribution (marginal)
#       InverseGammaDistribution (marginal)
#   2. (out):
#       Message{DeltaDistribution}
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
#       GaussianNode(name="my_node", V=1.0)
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
#
# GaussianNode with fixed mean and variance:
#
#       out
#   [N]----->
#
#   out = Message(GaussianDistribution(m=mean, V=variance))
#
#   Example:
#       GaussianNode(name="my_node", m=1.0, V=1.0)
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

    function GaussianNode(; name=unnamedStr(), form::ASCIIString="moment", m::Union(Float64,Vector{Float64},Nothing)=nothing, V::Union(Float64,Matrix{Float64},Nothing)=nothing)
        if m!=nothing && V!=nothing
            total_interfaces = 1
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

isDeterministic(::GaussianNode) = false


############################################
# Standard update functions
############################################

function sumProduct!{T1<:Any, T2<:Any}(node::GaussianNode,
                                       outbound_interface_id::Int,
                                       ::Nothing,
                                       msg_var_prec::Message{DeltaDistribution{T1}},
                                       msg_out::Message{DeltaDistribution{T2}})
    # Rules from Korl table 5.2 by symmetry
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if isdefined(node, :variance)
        # Backward over mean
        #
        #   N       Dlt            
        #  ---->[N]<----
        #  <--   |
        #        | Dlt
        #        v

        (length(msg_out.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_out.payload.m[1]]
        dist_out.xi = nothing
        dist_out.V = reshape([msg_var_prec.payload.m[1]], 1, 1)
        dist_out.W = nothing
    elseif isdefined(node, :precision)
        # Backward over mean
        #
        #   N       Dlt            
        #  ---->[N]<----
        #  <--   |
        #        | Dlt
        #        v

        (length(msg_out.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_out.payload.m[1]]
        dist_out.xi = nothing
        dist_out.V = nothing
        dist_out.W = reshape([msg_var_prec.payload.m[1]], 1, 1)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!{T1<:Any, T2<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            msg_mean::Message{DeltaDistribution{T1}},
                            ::Nothing,
                            msg_out::Message{DeltaDistribution{T2}})
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if isdefined(node, :variance)
        # Backward over variance
        #
        #   Dlt      IG            
        #  ---->[N]<----
        #        |   -->
        #     Dlt|  
        #        v

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload
        y = msg_out.payload.m[1]
        length(msg_mean.payload.m) == 1 || error("Update only defined for univariate distributions")
        m = msg_mean.payload.m[1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2
    elseif isdefined(node, :precision)
        # Backward over precision
        #
        #   Dlt      Gam            
        #  ---->[N]<----
        #        |   -->
        #     Dlt|  
        #        v

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload
        y = msg_out.payload.m[1]
        length(msg_mean.payload.m) == 1 || error("Update only defined for univariate distributions")
        m = msg_mean.payload.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y-m)^2
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!{T<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            msg_out::Message{DeltaDistribution{T}})
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if isdefined(node, :variance)
        # Backward over variance with fixed mean
        #
        #  Ig        Dlt               
        #  ---->[N]---->
        #   <--  

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload
        y = msg_out.payload.m[1]
        length(node.m) == 1 || error("Update only defined for univariate distributions")
        m = node.m[1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2
    elseif isdefined(node, :precision)
        # Backward over precision with fixed mean
        #
        #  Gam       Dlt               
        #  ---->[N]---->
        #   <--  

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload
        y = msg_out.payload.m[1]
        length(node.m) == 1 || error("Update only defined for univariate distributions")
        m = node.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y-m)^2
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!{T1<:Any, T2<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            msg_mean::Message{DeltaDistribution{T1}},
                            msg_var_prec::Message{DeltaDistribution{T2}},
                            ::Nothing)
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if isdefined(node, :variance)
        # Forward over out
        #
        #   Dlt     Dlt            
        #  ---->[N]<----
        #        |
        #      | |  
        #      v v

        (length(msg_mean.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_mean.payload.m[1]]
        dist_out.xi = nothing
        dist_out.V = reshape([msg_var_prec.payload.m[1]], 1, 1)
        dist_out.W = nothing
    elseif isdefined(node, :precision)
        # Forward over out
        #
        #   Dlt     Dlt            
        #  ---->[N]<----
        #        |
        #      | |  
        #      v v

        (length(msg_mean.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_mean.payload.m[1]]
        dist_out.xi = nothing
        dist_out.V = nothing
        dist_out.W = reshape([msg_var_prec.payload.m[1]], 1, 1)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!{T<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            msg_var_prec::Message{DeltaDistribution{T}},
                            msg_out::Nothing)
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if isdefined(node, :variance)
        # Forward over out with fixed mean
        #
        #  Dlt        N               
        #  ---->[N]---->
        #           -->  

        (length(node.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = deepcopy(node.m)
        dist_out.xi = nothing
        dist_out.V = reshape([msg_var_prec.payload.m[1]], 1, 1)
        dist_out.W = nothing
    elseif isdefined(node, :precision)
        # Forward over out with fixed mean
        #
        #  Dlt        N               
        #  ---->[N]---->
        #           -->  

        (length(node.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = deepcopy(node.m)
        dist_out.xi = nothing
        dist_out.V = nothing
        dist_out.W = reshape([msg_var_prec.payload.m[1]], 1, 1)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Any)

    # Forward over out with fixed mean and variance, effectively rendering the Gaussian node a terminal
    #
    #        N               
    #  [N]---->
    #      -->  

    return node.out.message = Message(GaussianDistribution(m=deepcopy(node.m), V=deepcopy(node.V)))
end

############################################
# Naive variational update functions
############################################

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            marg_variance::InverseGammaDistribution,
                            marg_out::GaussianDistribution)
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   <--     Q~IG            
    #  ---->[N]<----
    #    N   |   
    #        |Q~N  
    #        v

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.mean) # Mean estimation from variance and sample
        ensureMDefined!(marg_out)
        length(marg_out.m) == 1 || error("VMP for Gaussian node is only implemented for univariate distributions")
        a = marg_variance.a # gamma message
        b = marg_variance.b
        dist_out.m = deepcopy(marg_out.m)
        dist_out.V = reshape([((a+1)/b)], 1, 1)
        dist_out.xi = nothing
        dist_out.W = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
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

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

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

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            marg_out::GaussianDistribution)

    if isdefined(node, :precision) && is(node.interfaces[outbound_interface_id], node.precision)
        # Forward variational update function with fixed mean
        #
        #   Gam     Q~N                 
        #  ---->[N]---->
        #  <--

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload
        ensureMVParametrization!(marg_out)
        (length(node.m) == 1 && length(marg_out.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.a = 1.5
        dist_out.b = 0.5*(marg_out.m[1] - node.m[1])^2 + 0.5*marg_out.V[1,1]
    elseif isdefined(node, :mean) && is(node.interfaces[outbound_interface_id], node.mean)
        # Backward variational update function with fixed variance
        #
        #   <--     Q~N            
        #  ---->[N]---->
        #    N
        #

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload
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

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            marg_mean::GaussianDistribution,
                            ::Nothing,
                            marg_out::GaussianDistribution)
    if isdefined(node, :variance)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #   Q~N      IG            
        #  ---->[N]<----
        #        |   -->
        #    Q~N |  
        #        v

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload

        if is(node.interfaces[outbound_interface_id], node.variance) # Variance estimation from mean and sample
            ensureMVParametrization!(marg_out)
            ensureMVParametrization!(marg_mean)
            (length(marg_out.V) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
            (length(marg_mean.V) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
            mu_m = marg_mean.m[1]
            mu_y = marg_out.m[1]
            V_m = marg_mean.V[1,1]
            V_y = marg_out.V[1,1]
            dist_out.a = -0.5
            dist_out.b = 0.5*(mu_y-mu_m)^2+0.5*(V_m+V_y)
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    elseif isdefined(node, :precision)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #   Q~N      Gam            
        #  ---->[N]<----
        #       |   -->
        #   Q~N |  
        #       v

        dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload

        if is(node.interfaces[outbound_interface_id], node.precision) # Precision estimation from mean and sample
            ensureMVParametrization!(marg_out)
            ensureMVParametrization!(marg_mean)
            (length(marg_out.V) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
            (length(marg_mean.V) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
            mu_m = marg_mean.m[1]
            mu_y = marg_out.m[1]
            V_m = marg_mean.V[1,1]
            V_y = marg_out.V[1,1]
            dist_out.a = 1.5
            dist_out.b = 0.5*(mu_y-mu_m)^2+0.5*(V_m+V_y)
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    else
        error("Unknown update rule for $(typeof(node)) $(node.name). Only the moment and precision form are currently supported.")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
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

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

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

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            marg_prec::GammaDistribution,
                            ::Nothing)
    # Forward variational update function with fixed mean
    #
    #  Q~Gam      N                 
    #  ---->[N]---->
    #            -->

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.out)
        (length(node.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(node.m)
        dist_out.V = nothing
        dist_out.xi = nothing
        dist_out.W = reshape([marg_prec.a/marg_prec.b], 1, 1)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            marg_mean::GaussianDistribution,
                            ::Nothing)
    # Forward variational update function with fixed variance
    #
    #   Q~N     -->            
    #  ---->[N]---->
    #            N
    #

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

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

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
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

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

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


############################################
# Structured variational update functions
############################################

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
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

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], StudentsTDistribution).payload

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

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
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

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload

    if isdefined(node, :precision)
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

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
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

    dist_out = getOrCreateMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

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