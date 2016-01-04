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

type GaussianNode{node_form} <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    # Fixed parameters
    m::Vector{Float64}
    V::Matrix{Float64}

    function GaussianNode(id::Symbol, interfaces::Array{Interface,1}, i::Dict{Symbol,Interface})
        new(id, interfaces, i)
    end
end

function GaussianNode(; id=generateNodeId(GaussianNode), form::Symbol=:moment, m::Union{Float64,Vector{Float64}}=[NaN], V::Union{Float64,Matrix{Float64}}=reshape([NaN], 1, 1))
    if isValid(m) && isValid(V)
        total_interfaces = 1
    elseif isValid(m) || isValid(V)
        total_interfaces = 2
    else
        total_interfaces = 3
    end
    self = GaussianNode{Val{form}}(id, Array(Interface, total_interfaces), Dict{Symbol,Interface}())
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
        else
            error("Unrecognized form, $(form). Please use :moment, or :precision")
        end
        next_interface_index += 1
    end

    # Out interface
    self.i[:out] = self.interfaces[next_interface_index] = Interface(self)

    return self
end

isDeterministic(::GaussianNode) = false


############################################
# Sumproduct update functions
############################################

function sumProduct!(   node::GaussianNode{Val{:moment}},
                        outbound_interface_index::Type{Val{1}},
                        msg_mean::Any,
                        msg_var_prec::Message{DeltaDistribution{Float64}},
                        msg_out::Message{DeltaDistribution{Float64}},
                        outbound_dist::GaussianDistribution)

    # Rules from Korl table 5.2
    #
    #   N       Dlt
    #  ---->[N]<----
    #  <--   |
    #        | Dlt
    #        v

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var_prec.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:precision}},
                        outbound_interface_index::Type{Val{1}},
                        msg_mean::Any,
                        msg_var_prec::Message{DeltaDistribution{Float64}},
                        msg_out::Message{DeltaDistribution{Float64}},
                        outbound_dist::GaussianDistribution)

    # Rules from Korl table 5.2
    #
    #   N       Dlt
    #  ---->[N]<----
    #  <--   |
    #        | Dlt
    #        v

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_var_prec.payload.m

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:moment}},
                        outbound_interface_index::Type{Val{2}},
                        msg_mean::Message{DeltaDistribution{Float64}},
                        msg_var::Any,
                        msg_out::Message{DeltaDistribution{Float64}},
                        outbound_dist::InverseGammaDistribution)

    # Rules from Korl table 5.2
    #
    #   Dlt      IG
    #  ---->[N]<----
    #        |   -->
    #     Dlt|
    #        v

    y = msg_out.payload.m
    m = msg_mean.payload.m
    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:precision}},
                        outbound_interface_index::Type{Val{2}},
                        msg_mean::Message{DeltaDistribution{Float64}},
                        msg_prec::Any,
                        msg_out::Message{DeltaDistribution{Float64}},
                        outbound_dist::GammaDistribution)

    #   Dlt      Gam
    #  ---->[N]<----
    #        |   -->
    #     Dlt|
    #        v

    y = msg_out.payload.m
    m = msg_mean.payload.m
    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:moment}},
                        outbound_interface_index::Type{Val{3}},
                        msg_mean::Message{DeltaDistribution{Float64}},
                        msg_var_prec::Message{DeltaDistribution{Float64}},
                        msg_out::Any,
                        outbound_dist::GaussianDistribution)

    # Rules from Korl table 5.2
    #
    #   Dlt     Dlt
    #  ---->[N]<----
    #        |
    #      | | N
    #      v v

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var_prec.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:precision}},
                        outbound_interface_index::Type{Val{3}},
                        msg_mean::Message{DeltaDistribution{Float64}},
                        msg_var_prec::Message{DeltaDistribution{Float64}},
                        msg_out::Any,
                        outbound_dist::GaussianDistribution)

    # Rules from Korl table 5.2
    #
    #   Dlt     Dlt
    #  ---->[N]<----
    #        |
    #      | | N
    #      v v

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_var_prec.payload.m

    return outbound_dist
end


#############################################
# Sumproduct update functions with fixed mean
#############################################

function sumProduct!(   node::GaussianNode{Val{:moment}},
                        outbound_interface_index::Type{Val{1}},
                        msg_out::Any,
                        outbound_dist::GaussianDistribution)

    # Fixed mean and variance, effectively rendering the Gaussian node a terminal
    #
    #        N
    #  [N]---->
    #      -->

    outbound_dist.m = node.m[1]
    outbound_dist.V = node.V[1,1]
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:moment}},
                        outbound_interface_index::Type{Val{1}},
                        msg_var::Any,
                        msg_out::Message{DeltaDistribution{Float64}},
                        outbound_dist::InverseGammaDistribution)

    # Rules from Korl table 5.2
    # Fixed mean
    #
    #  Ig        Dlt
    #  ---->[N]---->
    #   <--

    y = msg_out.payload.m
    m = node.m[1]
    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:precision}},
                        outbound_interface_index::Type{Val{1}},
                        msg_prec::Any,
                        msg_out::Message{DeltaDistribution{Float64}},
                        outbound_dist::GammaDistribution)

    # Fixed mean
    #
    #  Gam       Dlt
    #  ---->[N]---->
    #   <--

    y = msg_out.payload.m
    m = node.m[1]
    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:moment}},
                        outbound_interface_index::Type{Val{2}},
                        msg_var::Message{DeltaDistribution{Float64}},
                        msg_out::Any,
                        outbound_dist::GaussianDistribution)

    # Rules from Korl table 5.2
    # Fixed mean
    #
    #  Dlt        N
    #  ---->[N]---->
    #           -->

    outbound_dist.m = node.m[1]
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProduct!(   node::GaussianNode{Val{:precision}},
                        outbound_interface_index::Type{Val{2}},
                        msg_prec::Message{DeltaDistribution{Float64}},
                        msg_out::Any,
                        outbound_dist::GaussianDistribution)

    # Rules from Korl table 5.2
    # Fixed mean
    #
    #  Dlt        N
    #  ---->[N]---->
    #           -->

    outbound_dist.m = node.m[1]
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_prec.payload.m

    return outbound_dist
end


############################################
# Naive variational update functions
############################################

function vmp!(  node::GaussianNode{Val{:moment}},
                outbound_interface_index::Type{Val{1}},
                marg_mean::Any,
                marg_variance::InverseGammaDistribution,
                marg_out::GaussianDistribution,
                outbound_dist::GaussianDistribution)

    # Derivation for the update rule can be found in the derivations notebook.
    #
    #   <--     Q~IG
    #  ---->[N]<----
    #    N   |
    #        |Q~N
    #        v

    (isProper(marg_variance) && isProper(marg_out)) || error("Improper input distributions are not supported")
    ensureParameters!(marg_out, (:m,))

    outbound_dist.m = marg_out.m
    outbound_dist.V = (marg_variance.a + 1)/marg_variance.b
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{1}},
                marg_mean::Any,
                marg_prec::GammaDistribution,
                marg_y::GaussianDistribution,
                outbound_dist::GaussianDistribution)

    #   N       Q~Gam
    #  ---->[N]<----
    #   <--  |
    #    Q~N |
    #        v

    (isProper(marg_prec) && isProper(marg_y)) || error("Improper input distributions are not supported")
    ensureParameters!(marg_y, (:m,))

    outbound_dist.m = marg_y.m
    outbound_dist.W = marg_prec.a / marg_prec.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:moment}},
                outbound_interface_index::Type{Val{2}},
                marg_mean::GaussianDistribution,
                marg_var::Any,
                marg_out::GaussianDistribution,
                outbound_dist::InverseGammaDistribution)

    #   Q~N      IG
    #  ---->[N]<----
    #        |   -->
    #    Q~N |
    #        v

    ensureParameters!(marg_out, (:m, :V))
    ensureParameters!(marg_mean, (:m, :V))

    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(marg_out.m - marg_mean.m)^2 + 0.5*(marg_mean.V + marg_out.V)

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{2}},
                marg_mean::GaussianDistribution,
                marg_prec::Any,
                marg_out::GaussianDistribution,
                outbound_dist::GammaDistribution)

    #   Q~N      Gam
    #  ---->[N]<----
    #        |   -->
    #    Q~N |
    #        v

    ensureParameters!(marg_out, (:m, :V))
    ensureParameters!(marg_mean, (:m, :V))

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(marg_out.m - marg_mean.m)^2 + 0.5*(marg_mean.V + marg_out.V)

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:moment}},
                outbound_interface_index::Type{Val{3}},
                marg_mean::GaussianDistribution,
                marg_var::InverseGammaDistribution,
                marg_out::Any,
                outbound_dist::GaussianDistribution)

    #   Q~N     Q~IG
    #  ---->[N]<----
    #        |
    #      N | |
    #        v v

    (isProper(marg_mean) && isProper(marg_var)) || error("Improper input distributions are not supported")
    ensureParameters!(marg_mean, (:m,))

    outbound_dist.m = marg_mean.m
    outbound_dist.V = marg_var.b / (marg_var.a-1.0)
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{3}},
                marg_mean::GaussianDistribution,
                marg_prec::GammaDistribution,
                marg_out::Any,
                outbound_dist::GaussianDistribution)

    #   Q~N     Q~Gam
    #  ---->[N]<----
    #        |
    #      N | |
    #        v v

    (isProper(marg_mean) && isProper(marg_prec)) || error("Improper input distributions are not supported")
    ensureParameters!(marg_mean, (:m,))

    outbound_dist.m = marg_mean.m
    outbound_dist.W = marg_prec.a / marg_prec.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end


##########################################################
# Naive variational update functions with fixed interfaces
##########################################################

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{1}},
                marg_prec::Any,
                marg_out::GaussianDistribution,
                outbound_dist::GammaDistribution)

    # Fixed mean
    #
    #   Gam     Q~N
    #  ---->[N]---->
    #  <--

    isProper(marg_out) || error("Improper input distributions are not supported")
    ensureParameters!(marg_out, (:m, :V))

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(marg_out.m - node.m[1])^2 + 0.5*marg_out.V

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:moment}},
                outbound_interface_index::Type{Val{1}},
                marg_mean::Any,
                marg_out::GaussianDistribution,
                outbound_dist::GaussianDistribution)

    # Fixed variance
    #
    #   <--     Q~N
    #  ---->[N]---->
    #    N

    isProper(marg_out) || error("Improper input distributions are not supported")
    ensureParameters!(marg_out, (:m,))

    outbound_dist.m = marg_out.m
    outbound_dist.V = node.V[1,1]
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{2}},
                marg_prec::GammaDistribution,
                marg_out::Any,
                outbound_dist::GaussianDistribution)

    # Fixed mean
    #
    #  Q~Gam      N
    #  ---->[N]---->
    #            -->

    isProper(marg_prec) || error("Improper input distributions are not supported")

    outbound_dist.m = node.m[1]
    outbound_dist.V = NaN
    outbound_dist.xi = NaN
    outbound_dist.W = marg_prec.a / marg_prec.b

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:moment}},
                outbound_interface_index::Type{Val{2}},
                marg_mean::GaussianDistribution,
                marg_out::Any,
                outbound_dist::GaussianDistribution)

    # Fixed variance
    #
    #   Q~N     -->
    #  ---->[N]---->
    #            N
    
    isProper(marg_mean) || error("Improper input distributions are not supported")
    ensureParameters!(marg_mean, (:m,))

    outbound_dist.m = marg_mean.m
    outbound_dist.V = node.V[1,1]
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end


############################################
# Structured variational update functions
############################################

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{1}},
                msg_mean::Any,
                msg_prec::Message{GammaDistribution},
                marg_out::GaussianDistribution,
                outbound_dist::StudentsTDistribution)

    #   <--     Gam
    #  ---->[N]<----
    #   St   |
    # - - - - - - - -
    #        | Q~N
    #        v

    (isProper(msg_prec.payload) && isProper(marg_out)) || error("Improper input distributions are not supported")
    ensureParameters!(marg_out, (:m, :W))

    outbound_dist.m = marg_out.m
    outbound_dist.lambda = (2.0*marg_out.W*msg_prec.payload.a)/(2.0*marg_out.W*msg_prec.payload.b + 1)
    outbound_dist.nu = 2.0*msg_prec.payload.a

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{2}},
                msg_mean::Message{GaussianDistribution},
                msg_prec::Any,
                marg_out::GaussianDistribution,
                outbound_dist::GammaDistribution)

    #    N      Gam
    #  ---->[N]<----
    #        |   -->
    # - - - - - - - -
    #        | Q~N
    #        v

    (isProper(msg_mean.payload) && isProper(marg_out)) || error("Improper input distributions are not supported")
    ensureParameters!(marg_out, (:m, :W))
    ensureParameters!(msg_mean.payload, (:m,))

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*((1.0/marg_out.W) + (msg_mean.payload.m - marg_out.m)^2)

    return outbound_dist
end

function vmp!(  node::GaussianNode{Val{:precision}},
                outbound_interface_index::Type{Val{3}},
                marg::NormalGammaDistribution,
                ::NormalGammaDistribution, # Same distribution as marg
                marg_out::Any,
                outbound_dist::GaussianDistribution)

    #      Q~NGam
    #  ---->[N]<----
    #        |
    # - - - - - - - -
    #      | | N
    #      v v

    isProper(marg) || error("Improper input distributions are not supported")

    outbound_dist.m = marg.m
    outbound_dist.W = marg.a / marg.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end
