############################################
# GaussianNode
############################################
# Description:
#   Node converting an input mean and precision
#   to a univariate Gaussian distribution:
#   out = N(mean, variance) or
#   out = N(mean, precision⁻¹)
#   
#   The GaussianNode has two different roles for its second interface,
#   namely precision and variance. This is determined by the form argument,
#   either :moment or :precision.
#
#   Also the GaussianNode accepts a fixed mean and/or variance. In this case
#   the fixed interface(s) are not constructed and the interface indices shift
#   a position.
#
#           mean
#            |
#            v  out
#     ----->[N]----->
#    precision/
#    variance
#
# Interfaces:
#   1 i[:mean], 2 i[:variance] or i[:precision], 3 i[:out]
#
# Construction:
#   GaussianNode(form=:moment, id=:my_node, m=optional, V=optional)
#
############################################

export GaussianNode

type GaussianNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    # Fixed parameters
    m::Vector{Float64}
    V::Matrix{Float64}

    function GaussianNode(; id=generateNodeId(GaussianNode), form::Symbol=:moment, m::Union(Float64,Vector{Float64})=[NaN], V::Union(Float64,Matrix{Float64})=reshape([NaN], 1, 1))
        if isValid(m) && isValid(V)
            total_interfaces = 1
        elseif isValid(m) || isValid(V)
            total_interfaces = 2
        else
            total_interfaces = 3
        end
        self = new(id, Array(Interface, total_interfaces), Dict{Symbol,Interface}())
        !haskey(current_graph.n, id) ? current_graph.n[id] = self : error("Node id $(id) already present")
        next_interface_index = 1 # Counter keeping track of constructed interfaces

        if isValid(m)
            # GaussianNode with fixed mean
            self.m = (typeof(m)==Float64) ? [m] : deepcopy(m)
        else
            # GaussianNode with variable mean
            self.i[:mean] = self.interfaces[next_interface_index] = Interface(self) # Mean interface
            next_interface_index += 1
        end

        # Pick a form for the variance/precision
        if isValid(V)
            # GaussianNode with fixed variance
            self.V = (typeof(V)==Float64) ? reshape([V], 1, 1) : deepcopy(V)
        else
            # GaussianNode with variable variance
            self.interfaces[next_interface_index] = Interface(self)
            if form == :moment
                # Parameters m, V
                self.i[:variance] = self.interfaces[next_interface_index]
            elseif form == :precision
                # Parameters m, W
                self.i[:precision] = self.interfaces[next_interface_index]
            else
                error("Unrecognized form, $(form). Please use :moment, or :precision")
            end
            next_interface_index += 1
        end
        
        # Out interface
        self.i[:out] = self.interfaces[next_interface_index] = Interface(self)

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

    if haskey(node.i, :variance)
        # Backward over mean
        #
        #   N       Dlt            
        #  ---->[N]<----
        #  <--   |
        #        | Dlt
        #        v

        (length(msg_out.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_out.payload.m[1]]
        invalidate!(dist_out.xi)
        dist_out.V = reshape([msg_var_prec.payload.m[1]], 1, 1)
        invalidate!(dist_out.W)

        return (:gaussian_backward_mean_delta_variance,
                node.interfaces[outbound_interface_id].message)
    elseif haskey(node.i, :precision)
        # Backward over mean
        #
        #   N       Dlt            
        #  ---->[N]<----
        #  <--   |
        #        | Dlt
        #        v

        (length(msg_out.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_out.payload.m[1]]
        invalidate!(dist_out.xi)
        invalidate!(dist_out.V)
        dist_out.W = reshape([msg_var_prec.payload.m[1]], 1, 1)

        return (:gaussian_backward_mean_delta_precision,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!{T1<:Any, T2<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            msg_mean::Message{DeltaDistribution{T1}},
                            ::Nothing,
                            msg_out::Message{DeltaDistribution{T2}})
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if haskey(node.i, :variance)
        # Backward over variance
        #
        #   Dlt      IG            
        #  ---->[N]<----
        #        |   -->
        #     Dlt|  
        #        v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload
        y = msg_out.payload.m[1]
        length(msg_mean.payload.m) == 1 || error("Update only defined for univariate distributions")
        m = msg_mean.payload.m[1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_variance_delta,
                node.interfaces[outbound_interface_id].message)
    elseif haskey(node.i, :precision)
        # Backward over precision
        #
        #   Dlt      Gam            
        #  ---->[N]<----
        #        |   -->
        #     Dlt|  
        #        v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload
        y = msg_out.payload.m[1]
        length(msg_mean.payload.m) == 1 || error("Update only defined for univariate distributions")
        m = msg_mean.payload.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_precision_delta,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!{T<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            msg_out::Message{DeltaDistribution{T}})
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if haskey(node.i, :variance)
        # Backward over variance with fixed mean
        #
        #  Ig        Dlt               
        #  ---->[N]---->
        #   <--  

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload
        y = msg_out.payload.m[1]
        length(node.m) == 1 || error("Update only defined for univariate distributions")
        m = node.m[1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_variance_delta,
                node.interfaces[outbound_interface_id].message)
    elseif haskey(node.i, :precision)
        # Backward over precision with fixed mean
        #
        #  Gam       Dlt               
        #  ---->[N]---->
        #   <--  

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload
        y = msg_out.payload.m[1]
        length(node.m) == 1 || error("Update only defined for univariate distributions")
        m = node.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_precision_delta,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!{T1<:Any, T2<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            msg_mean::Message{DeltaDistribution{T1}},
                            msg_var_prec::Message{DeltaDistribution{T2}},
                            ::Nothing)
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if haskey(node.i, :variance)
        # Forward over out
        #
        #   Dlt     Dlt            
        #  ---->[N]<----
        #        |
        #      | |  
        #      v v

        (length(msg_mean.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_mean.payload.m[1]]
        invalidate!(dist_out.xi)
        dist_out.V = reshape([msg_var_prec.payload.m[1]], 1, 1)
        invalidate!(dist_out.W)

        return (:gaussian_forward_delta_variance,
                node.interfaces[outbound_interface_id].message)
    elseif haskey(node.i, :precision)
        # Forward over out
        #
        #   Dlt     Dlt            
        #  ---->[N]<----
        #        |
        #      | |  
        #      v v

        (length(msg_mean.payload.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = [msg_mean.payload.m[1]]
        invalidate!(dist_out.xi)
        invalidate!(dist_out.V)
        dist_out.W = reshape([msg_var_prec.payload.m[1]], 1, 1)

        return (:gaussian_forward_delta_precision,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!{T<:Any}(node::GaussianNode,
                            outbound_interface_id::Int,
                            msg_var_prec::Message{DeltaDistribution{T}},
                            msg_out::Nothing)
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if haskey(node.i, :variance)
        # Forward over out with fixed mean
        #
        #  Dlt        N               
        #  ---->[N]---->
        #           -->  

        (length(node.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = deepcopy(node.m)
        invalidate!(dist_out.xi)
        dist_out.V = reshape([msg_var_prec.payload.m[1]], 1, 1)
        invalidate!(dist_out.W)

        return (:gaussian_forward_delta_variance,
                node.interfaces[outbound_interface_id].message)
    elseif haskey(node.i, :precision)
        # Forward over out with fixed mean
        #
        #  Dlt        N               
        #  ---->[N]---->
        #           -->  

        (length(node.m) == 1) || error("GaussianNode with fixed mean update only implemented for unvariate distributions")
        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

        dist_out.m = deepcopy(node.m)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.V)
        dist_out.W = reshape([msg_var_prec.payload.m[1]], 1, 1)

        return (:gaussian_forward_delta_precision,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Any)

    # Forward over out with fixed mean and variance, effectively rendering the Gaussian node a terminal
    #
    #        N               
    #  [N]---->
    #      -->  
    node.i[:out].message = Message(GaussianDistribution(m=deepcopy(node.m), V=deepcopy(node.V)))
    return (:gaussian_forward_fixed_mean_variance,
            node.interfaces[outbound_interface_id].message)
end

############################################
# Naive variational update functions
############################################

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:mean]) # Mean estimation from variance and sample
        ensureMDefined!(marg_out)
        length(marg_out.m) == 1 || error("VMP for Gaussian node is only implemented for univariate distributions")
        a = marg_variance.a # gamma message
        b = marg_variance.b
        dist_out.m = deepcopy(marg_out.m)
        dist_out.V = reshape([((a+1)/b)], 1, 1)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.W)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end

    return (:gaussian_backward_mean_gaussian_gamma,
            node.interfaces[outbound_interface_id].message)
end

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:mean])
        ensureMDefined!(marg_y)
        (length(marg_y.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_y.m)
        dist_out.W = reshape([marg_prec.a/marg_prec.b], 1, 1)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.V)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end

    return (:gaussian_backward_mean_gaussian_inverse_gamma,
            node.interfaces[outbound_interface_id].message)
end

function vmp!(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            marg_out::GaussianDistribution)

    if haskey(node.i, :precision) && is(node.interfaces[outbound_interface_id], node.i[:precision])
        # Forward variational update function with fixed mean
        #
        #   Gam     Q~N                 
        #  ---->[N]---->
        #  <--

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload
        ensureMVParametrization!(marg_out)
        (length(node.m) == 1 && length(marg_out.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.a = 1.5
        dist_out.b = 0.5*(marg_out.m[1] - node.m[1])^2 + 0.5*marg_out.V[1,1]

        return (:gaussian_backward_precision_gaussian_delta,
                node.interfaces[outbound_interface_id].message)
    elseif haskey(node.i, :mean) && is(node.interfaces[outbound_interface_id], node.i[:mean])
        # Backward variational update function with fixed variance
        #
        #   <--     Q~N            
        #  ---->[N]---->
        #    N
        #

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload
        ensureMDefined!(marg_out)
        dist_out.m = deepcopy(marg_out.m)
        dist_out.V = deepcopy(node.V)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.W)

        return (:gaussian_backward_mean_gaussian_delta,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(node::GaussianNode,
                            outbound_interface_id::Int,
                            marg_mean::GaussianDistribution,
                            ::Nothing,
                            marg_out::GaussianDistribution)
    if haskey(node.i, :variance)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #   Q~N      IG            
        #  ---->[N]<----
        #        |   -->
        #    Q~N |  
        #        v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], InverseGammaDistribution).payload

        if is(node.interfaces[outbound_interface_id], node.i[:variance]) # Variance estimation from mean and sample
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

            return (:gaussian_backward_variance_gaussian,
                    node.interfaces[outbound_interface_id].message)
        else
            error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
        end
    elseif haskey(node.i, :precision)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #   Q~N      Gam            
        #  ---->[N]<----
        #       |   -->
        #   Q~N |  
        #       v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload

        if is(node.interfaces[outbound_interface_id], node.i[:precision]) # Precision estimation from mean and sample
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

            return (:gaussian_backward_precision_gaussian,
                    node.interfaces[outbound_interface_id].message)
        else
            error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
        end
    else
        error("Unknown update rule for $(typeof(node)) $(node.id). Only the moment and precision form are currently supported.")
    end
end

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:out])
        ensureMDefined!(marg_mean)
        (length(marg_mean.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_mean.m)
        dist_out.V = reshape([marg_var.b/(marg_var.a-1.0)], 1, 1)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.W)

        return (:gaussian_forward_gaussian_inverse_gamma,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(node::GaussianNode,
                            outbound_interface_id::Int,
                            marg_prec::GammaDistribution,
                            ::Nothing)
    # Forward variational update function with fixed mean
    #
    #  Q~Gam      N                 
    #  ---->[N]---->
    #            -->

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:out])
        (length(node.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(node.m)
        invalidate!(dist_out.V)
        invalidate!(dist_out.xi)
        dist_out.W = reshape([marg_prec.a/marg_prec.b], 1, 1)

        return (:gaussian_forward_gamma,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(node::GaussianNode,
                            outbound_interface_id::Int,
                            marg_mean::GaussianDistribution,
                            ::Nothing)
    # Forward variational update function with fixed variance
    #
    #   Q~N     -->            
    #  ---->[N]---->
    #            N
    #

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:out])
        ensureMDefined!(marg_mean)
        dist_out.m = deepcopy(marg_mean.m)
        dist_out.V = deepcopy(node.V)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.W)

        return (:gaussian_forward_gaussian,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:out])
        ensureMDefined!(marg_mean)
        (length(marg_mean.m) == 1) || error("VMP for Gaussian node is only implemented for univariate distributions")
        dist_out.m = deepcopy(marg_mean.m)
        dist_out.W = reshape([marg_prec.a/marg_prec.b], 1, 1)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.V)

        return (:gaussian_forward_gaussian_gamma,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end


############################################
# Structured variational update functions
############################################

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], StudentsTDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:mean])
        ensureMWParametrization!(marg_out)
        length(marg_out.m) == 1 || error("SVMP update on GaussianNode only defined for univariate distributions")
        m_y = marg_out.m[1]
        prec_y = marg_out.W[1, 1]
        a = msg_prec.payload.a
        b = msg_prec.payload.b
        dist_out.m = [m_y]
        dist_out.W = reshape([(2.0*prec_y*a)/(2.0*prec_y*b + 1)], 1, 1)
        dist_out.nu = 2.0*a

        return (:gaussian_backward_mean_structured,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GammaDistribution).payload

    if haskey(node.i, :precision)
        ensureMWParametrization!(marg_out)
        length(marg_out.m) == 1 || error("SVMP update on GaussianNode only defined for univariate distributions")
        m_y = marg_out.m[1]
        prec_y = marg_out.W[1, 1]
        ensureMDefined!(msg_mean.payload)
        m_m = msg_mean.payload.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*((1.0/prec_y) + (m_m - m_y)^2)

        return (:gaussian_backward_precision_structured,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(node::GaussianNode,
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

    dist_out = ensureMessage!(node.interfaces[outbound_interface_id], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_id], node.i[:out])
        dist_out.m = [marg.m]
        dist_out.W = reshape([marg.a/marg.b], 1, 1)
        invalidate!(dist_out.xi)
        invalidate!(dist_out.V)

        return (:gaussian_forward_structured,
                node.interfaces[outbound_interface_id].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end