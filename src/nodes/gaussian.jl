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
#    log-variance/
#    variance
#
# Interfaces:
#   1 i[:mean], 2 i[:precision] / i[:log_variance] / i[:variance], 3 i[:out]
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

    function GaussianNode(; id=generateNodeId(GaussianNode), form::Symbol=:moment, m::Union{Float64,Vector{Float64}}=[NaN], V::Union{Float64,Matrix{Float64}}=reshape([NaN], 1, 1))
        if isValid(m) && isValid(V)
            total_interfaces = 1
        elseif isValid(m) || isValid(V)
            total_interfaces = 2
        else
            total_interfaces = 3
        end
        self = new(id, Array(Interface, total_interfaces), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

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
            elseif form == :log_variance
                # Parameters m, log(W)
                self.i[:log_variance] = self.interfaces[next_interface_index]
            else
                error("Unrecognized form, $(form). Please use :moment, :precision of :log_variance")
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

function sumProduct!(node::GaussianNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_var_prec::Message{DeltaDistribution{Float64}},
                     msg_out::Message{DeltaDistribution{Float64}})
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

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        dist_out.m = msg_out.payload.m
        dist_out.xi = NaN
        dist_out.V = msg_var_prec.payload.m
        dist_out.W = NaN

        return (:gaussian_backward_mean_delta_variance,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :precision)
        # Backward over mean
        #
        #   N       Dlt
        #  ---->[N]<----
        #  <--   |
        #        | Dlt
        #        v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        dist_out.m = msg_out.payload.m
        dist_out.xi = NaN
        dist_out.V = NaN
        dist_out.W = msg_var_prec.payload.m

        return (:gaussian_backward_mean_delta_precision,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!(node::GaussianNode,
                     outbound_interface_index::Int,
                     msg_mean::Message{DeltaDistribution{Float64}},
                     ::Void,
                     msg_out::Message{DeltaDistribution{Float64}})
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

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], InverseGammaDistribution).payload
        y = msg_out.payload.m
        m = msg_mean.payload.m
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_variance_delta,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :precision)
        # Backward over precision
        #
        #   Dlt      Gam
        #  ---->[N]<----
        #        |   -->
        #     Dlt|
        #        v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GammaDistribution).payload
        y = msg_out.payload.m
        m = msg_mean.payload.m
        dist_out.a = 1.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_precision_delta,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!(node::GaussianNode,
                     outbound_interface_index::Int,
                     ::Void,
                     msg_out::Message{DeltaDistribution{Float64}})
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if haskey(node.i, :variance)
        # Backward over variance with fixed mean
        #
        #  Ig        Dlt
        #  ---->[N]---->
        #   <--

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], InverseGammaDistribution).payload
        y = msg_out.payload.m
        m = node.m[1]
        dist_out.a = -0.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_variance_delta,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :precision)
        # Backward over precision with fixed mean
        #
        #  Gam       Dlt
        #  ---->[N]---->
        #   <--

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GammaDistribution).payload
        y = msg_out.payload.m
        m = node.m[1]
        dist_out.a = 1.5
        dist_out.b = 0.5*(y-m)^2

        return (:gaussian_backward_precision_delta,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!(   node::GaussianNode,
                        outbound_interface_index::Int,
                        msg_mean::Message{DeltaDistribution{Float64}},
                        msg_var_prec::Message{DeltaDistribution{Float64}},
                        ::Void)
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

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        dist_out.m = msg_mean.payload.m
        dist_out.xi = NaN
        dist_out.V = msg_var_prec.payload.m
        dist_out.W = NaN

        return (:gaussian_forward_delta_variance,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :precision)
        # Forward over out
        #
        #   Dlt     Dlt
        #  ---->[N]<----
        #        |
        #      | |
        #      v v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        dist_out.m = msg_mean.payload.m
        dist_out.xi = NaN
        dist_out.V = NaN
        dist_out.W = msg_var_prec.payload.m

        return (:gaussian_forward_delta_precision,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!(   node::GaussianNode,
                        outbound_interface_index::Int,
                        msg_var_prec::Message{DeltaDistribution{Float64}},
                        msg_out::Void)
    # Rules from Korl table 5.2
    # Note that in the way Korl wrote this it is an approximation; the actual result would be a student's t.
    # Here we assume the variance to be a point estimate.

    if haskey(node.i, :variance)
        # Forward over out with fixed mean
        #
        #  Dlt        N
        #  ---->[N]---->
        #           -->

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        dist_out.m = node.m[1]
        dist_out.xi = NaN
        dist_out.V = msg_var_prec.payload.m
        dist_out.W = NaN

        return (:gaussian_forward_delta_variance,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :precision)
        # Forward over out with fixed mean
        #
        #  Dlt        N
        #  ---->[N]---->
        #           -->

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        dist_out.m = node.m[1]
        dist_out.xi = NaN
        dist_out.V = NaN
        dist_out.W = msg_var_prec.payload.m

        return (:gaussian_forward_delta_precision,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function sumProduct!(   node::GaussianNode,
                        outbound_interface_index::Int,
                        ::Any)

    # Forward over out with fixed mean and variance, effectively rendering the Gaussian node a terminal
    #
    #        N
    #  [N]---->
    #      -->
    node.i[:out].message = Message(GaussianDistribution(m=node.m[1], V=node.V[1,1]))
    return (:gaussian_forward_fixed_mean_variance,
            node.interfaces[outbound_interface_index].message)
end

############################################
# Naive variational update functions
############################################

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                ::Void,
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
    # (isProper(marg_variance) && isProper(marg_out)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:mean]) # Mean estimation from variance and sample
        ensureParameters!(marg_out, (:m,))
        a = marg_variance.a # gamma message
        b = marg_variance.b
        dist_out.m = marg_out.m
        dist_out.V = (a+1)/b
        dist_out.xi = NaN
        dist_out.W = NaN
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end

    return (:gaussian_backward_mean_gaussian_inverse_gamma,
            node.interfaces[outbound_interface_index].message)
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                ::Void,
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
    # (isProper(marg_prec) && isProper(marg_y)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:mean])
        ensureParameters!(marg_y, (:m,))
        dist_out.m = marg_y.m
        dist_out.W = marg_prec.a / marg_prec.b
        dist_out.xi = NaN
        dist_out.V = NaN
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end

    return (:gaussian_backward_mean_gaussian_gamma,
            node.interfaces[outbound_interface_index].message)
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                ::Void,
                marg_log_prec::GaussianDistribution,
                marg_y::GaussianDistribution)
    # Backward over mean, Gaussian variance input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # By symmetry the rules are the same as the forward over out.
    #
    #   N       Q~N (log-variance)
    #  ---->[N]<----
    #   <--  |
    #    Q~N |
    #        v
    # (isProper(marg_log_prec) && isProper(marg_y)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:mean])
        ensureParameters!(marg_y, (:m,))
        ensureParameters!(marg_log_prec, (:m, :V))
        dist_out.m = marg_y.m
        V = exp(marg_log_prec.m + 0.5*marg_log_prec.V)
        dist_out.V = (V > huge ? huge : V)
        dist_out.xi = NaN
        dist_out.W = NaN
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end

    return (:gaussian_backward_mean_gaussian_log_variance,
            node.interfaces[outbound_interface_index].message)
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                ::Void,
                marg_out::GaussianDistribution)
    # isProper(marg_out) || error("Improper input distributions are not supported")
    if haskey(node.i, :precision) && is(node.interfaces[outbound_interface_index], node.i[:precision])
        # Backward variational update function with fixed mean
        #
        #   Gam     Q~N
        #  ---->[N]---->
        #  <--

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GammaDistribution).payload
        ensureParameters!(marg_out, (:m, :V))
        dist_out.a = 1.5
        dist_out.b = 0.5*(marg_out.m - node.m[1])^2 + 0.5*marg_out.V

        return (:gaussian_backward_precision_gaussian_delta,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :variance) && is(node.interfaces[outbound_interface_index], node.i[:variance])
        # Backward variational update function with fixed mean
        #
        #   Ig      Q~N
        #  ---->[N]---->
        #  <--

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], InverseGammaDistribution).payload
        ensureParameters!(marg_out, (:m, :V))
        dist_out.a = -0.5
        dist_out.b = 0.5*(marg_out.m - node.m[1])^2 + 0.5*marg_out.V

        return (:gaussian_backward_variance_gaussian_delta,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :log_variance) && is(node.interfaces[outbound_interface_index], node.i[:log_variance])
        # Backward variational update function with fixed mean
        #
        #    N      Q~N
        #  ---->[N]---->
        #  <--

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
        ensureParameters!(marg_out, (:m, :V))
        dist_out.m = log(node.m[1]^2 - 2.0*node.m[1]*marg_out.m + marg_out.m^2 + marg_out.V)
        dist_out.V = 2.0
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_backward_log_variance_gaussian_delta,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :mean) && is(node.interfaces[outbound_interface_index], node.i[:mean])
        # Backward variational update function with fixed variance
        #
        #   <--     Q~N
        #  ---->[N]---->
        #    N
        #

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
        ensureParameters!(marg_out, (:m,))
        dist_out.m = marg_out.m
        dist_out.V = node.V[1,1]
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_backward_mean_gaussian_delta,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_mean::GaussianDistribution,
                ::Void,
                marg_out::GaussianDistribution)
    # (isProper(marg_mean) && isProper(marg_out)) || error("Improper input distributions are not supported")
    if haskey(node.i, :variance)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #   Q~N      IG
        #  ---->[N]<----
        #        |   -->
        #    Q~N |
        #        v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], InverseGammaDistribution).payload

        if is(node.interfaces[outbound_interface_index], node.i[:variance]) # Variance estimation from mean and sample
            ensureParameters!(marg_out, (:m, :V))
            ensureParameters!(marg_mean, (:m, :V))
            mu_m = marg_mean.m
            mu_y = marg_out.m
            V_m = marg_mean.V
            V_y = marg_out.V
            dist_out.a = -0.5
            dist_out.b = 0.5*(mu_y - mu_m)^2 + 0.5*(V_m + V_y)

            return (:gaussian_backward_variance_gaussian,
                    node.interfaces[outbound_interface_index].message)
        else
            error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
        end
    elseif haskey(node.i, :precision)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #  Q~N      Gam
        # ---->[N]<----
        #       |   -->
        #   Q~N |
        #       v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GammaDistribution).payload

        if is(node.interfaces[outbound_interface_index], node.i[:precision]) # Precision estimation from mean and sample
            ensureParameters!(marg_out, (:m, :V))
            ensureParameters!(marg_mean, (:m, :V))
            mu_m = marg_mean.m
            mu_y = marg_out.m
            V_m = marg_mean.V
            V_y = marg_out.V
            dist_out.a = 1.5
            dist_out.b = 0.5*(mu_y-mu_m)^2+0.5*(V_m+V_y)

            return (:gaussian_backward_precision_gaussian,
                    node.interfaces[outbound_interface_index].message)
        else
            error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
        end
    elseif haskey(node.i, :log_variance)
        # Variational update function takes the marginals as input (instead of the inbound messages)
        # Derivation for the update rule can be found in the derivations notebook.
        #
        #  Q~N      N (log-variance)
        # ---->[N]<----
        #       |   -->
        #   Q~N |
        #       v

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

        if is(node.interfaces[outbound_interface_index], node.i[:log_variance])
            ensureParameters!(marg_out, (:m, :V))
            ensureParameters!(marg_mean, (:m, :V))
            dist_out.m = log(marg_mean.m^2 + marg_mean.V - 2.0*marg_mean.m*marg_out.m + marg_out.m^2 + marg_out.V)
            dist_out.V = 2.0
            dist_out.xi = NaN
            dist_out.W = NaN

            return (:gaussian_backward_log_variance_gaussian,
                    node.interfaces[outbound_interface_index].message)
        else
            error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
        end
    else
        error("Unknown update rule for $(typeof(node)) $(node.id). Only the moment and precision form are currently supported.")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_mean::GaussianDistribution,
                marg_var::InverseGammaDistribution,
                ::Void)
    # Forward over out, InverseGamma input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N     Q~IG
    #  ---->[N]<----
    #        |
    #      N | |
    #        v v
    # (isProper(marg_mean) && isProper(marg_var)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:out])
        ensureParameters!(marg_mean, (:m,))
        dist_out.m = marg_mean.m
        dist_out.V = marg_var.b / (marg_var.a-1.0)
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_forward_gaussian_inverse_gamma,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_prec::GammaDistribution,
                ::Void)
    # Forward variational update function with fixed mean
    #
    #  Q~Gam      N
    #  ---->[N]---->
    #            -->
    # isProper(marg_prec) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:out])
        dist_out.m = node.m[1]
        dist_out.V = NaN
        dist_out.xi = NaN
        dist_out.W = marg_prec.a / marg_prec.b

        return (:gaussian_forward_gamma,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_var::InverseGammaDistribution,
                ::Void)
    # Forward variational update function with fixed mean
    #
    #  Q~Ig      N
    #  ---->[N]---->
    #            -->
    # isProper(marg_var) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:out])
        dist_out.m = node.m[1]
        dist_out.V = marg_var.b/(marg_var.a - 1.0)
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_forward_inverse_gamma,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_in::GaussianDistribution,
                ::Void)

    # isProper(marg_in) || error("Improper input distributions are not supported")

    if haskey(node.i, :log_variance) && is(node.interfaces[outbound_interface_index], node.i[:out])
        # Forward variational update function with fixed mean
        #
        #   Q~N     -->
        #  ---->[N]---->
        #            N
        #

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
        ensureParameters!(marg_in, (:m,))
        dist_out.m = node.m[1]
        V = exp(marg_in.m + 0.5*marg_in.V)
        dist_out.V = (V > huge ? huge : V)
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_forward_log_variance_gaussian,
                node.interfaces[outbound_interface_index].message)
    elseif haskey(node.i, :mean) && is(node.interfaces[outbound_interface_index], node.i[:out])
        # Forward variational update function with fixed variance
        #
        #   Q~N     -->
        #  ---->[N]---->
        #            N
        #

        dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload
        ensureParameters!(marg_in, (:m,))
        dist_out.m = marg_in.m
        dist_out.V = node.V[1,1]
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_forward_gaussian,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_mean::GaussianDistribution,
                marg_prec::GammaDistribution,
                ::Void)
    # Forward over out, Gamma input
    # Variational update function takes the marginals as input (instead of the inbound messages)
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N     Q~Gam
    #  ---->[N]<----
    #        |
    #      N | |
    #        v v
    # (isProper(marg_mean) && isProper(marg_prec)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:out])
        ensureParameters!(marg_mean, (:m,))
        dist_out.m = marg_mean.m
        dist_out.W = marg_prec.a / marg_prec.b
        dist_out.xi = NaN
        dist_out.V = NaN

        return (:gaussian_forward_gaussian_gamma,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg_mean::GaussianDistribution,
                marg_log_prec::GaussianDistribution,
                ::Void)
    # Forward over out, Gaussian input
    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   Q~N     Q~N (log-variance)
    #  ---->[N]<----
    #        |
    #      N | |
    #        v v
    # (isProper(marg_mean) && isProper(marg_log_prec)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:out])
        ensureParameters!(marg_mean, (:m,))
        ensureParameters!(marg_log_prec, (:m,:V))
        dist_out.m = marg_mean.m
        V = exp(marg_log_prec.m + 0.5*marg_log_prec.V)
        dist_out.V = (V > huge ? huge : V)
        dist_out.xi = NaN
        dist_out.W = NaN

        return (:gaussian_forward_gaussian_log_variance,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

############################################
# Structured variational update functions
############################################

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                ::Void,
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
    # (isProper(msg_prec.payload) && isProper(marg_out)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], StudentsTDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:mean])
        ensureParameters!(marg_out, (:m, :W))
        m_y = marg_out.m
        prec_y = marg_out.W
        a = msg_prec.payload.a
        b = msg_prec.payload.b
        dist_out.m = m_y
        dist_out.lambda = (2.0*prec_y*a) / (2.0*prec_y*b + 1)
        dist_out.nu = 2.0*a

        return (:gaussian_backward_mean_structured,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                msg_mean::Message{GaussianDistribution},
                ::Void,
                marg_out::GaussianDistribution)
    # Backward message over the precision
    #
    #    N      Gam
    #  ---->[N]<----
    #        |   -->
    # - - - - - - - -
    #        | Q~N
    #        v
    # (isProper(msg_mean.payload) && isProper(marg_out)) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GammaDistribution).payload

    if haskey(node.i, :precision)
        ensureParameters!(marg_out, (:m, :W))
        m_y = marg_out.m
        prec_y = marg_out.W
        ensureParameters!(msg_mean.payload, (:m,))
        m_m = msg_mean.payload.m
        dist_out.a = 1.5
        dist_out.b = 0.5*((1.0/prec_y) + (m_m - m_y)^2)

        return (:gaussian_backward_precision_structured,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end

function vmp!(  node::GaussianNode,
                outbound_interface_index::Int,
                marg::NormalGammaDistribution,
                ::NormalGammaDistribution, # Same distribution as marg
                ::Void)
    # Forward message over out.
    # This update function has two argments instead of three because it uses the node's joint marginal over the mean and precision.
    #
    #      Q~NGam
    #  ---->[N]<----
    #        |
    # - - - - - - - -
    #      | | N
    #      v v
    # isProper(marg) || error("Improper input distributions are not supported")
    dist_out = ensureMessage!(node.interfaces[outbound_interface_index], GaussianDistribution).payload

    if is(node.interfaces[outbound_interface_index], node.i[:out])
        dist_out.m = marg.m
        dist_out.W = marg.a / marg.b
        dist_out.xi = NaN
        dist_out.V = NaN

        return (:gaussian_forward_structured,
                node.interfaces[outbound_interface_index].message)
    else
        error("Undefined inbound-outbound message type combination for node $(node.id) of type $(typeof(node)).")
    end
end
