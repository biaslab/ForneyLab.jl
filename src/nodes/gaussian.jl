export GaussianNode

"""
Description:

    Gaussian factor node.

Interfaces:

    The GaussianNode can take multiple parametrization forms, specified
    through the form argument:

    - form == :variance     => f(out,m,V)  = N(out|m,V)        (default)
    - form == :precision    => f(out,m,W)  = N(out|m,inv(W))
    - form == :log_variance => f(out,m,lV) = N(out|m,exp(lV))

    Each parameter can be fixed, in which case the corresponding interface does
    not exist (and the ids of subsequent interfaces are lowered).
    Examples:

    No fixed parameters:            |     Fixed mean:
    --------------------            |     -----------
                                    |
           1. :mean                 |
             |                      |
             v                      |
    2. ----->[N]-----> 3. :out      |     1. ----->[N]-----> 2. :out
     :variance                      |       :precision
                                    |
    GaussianNode(form=:variance)    |     GaussianNode(m=0.0, form=:precision)
                                    |
    ________________________________|______________________________________________________
                                    |
    Fixed variance:                 |     Fixed mean and variance (constant distribution):
    ---------------                 |     ------------------------------------------------
                                    |
                                    |
    1. ----->[N]-----> 2. :out      |     [N]-----> 1. :out
     :mean                          |
                                    |
    GaussianNode(V=1.0)             |     GaussianNode(m=0.0, V=1.0)
                                    |

    The type is parametrized: GaussianNode{mean_type,uncertainty_type}, where
       - mean_type ∈ {:fixed_mean, :mean}
       - uncertainty_type ∈ {:fixed_variance, :variance, :precision, :log_variance}

Construction:

    GaussianNode(form=:variance, id=:my_node)           # N(out|m,V)        (default)
    GaussianNode(form=:precision, m=0.0, id=:my_node)   # N(out|0.0,inv(W)) (fixed mean)
    GaussianNode(V=1.0, id=:my_node)                    # N(out|m,1.0)      (fixed variance)
    GaussianNode(m=0.0, V=1.0, id=:my_node)             # N(out|0.0,1.0)    (fixed distribution)
"""
type GaussianNode{mean_type,uncertainty_type} <: Node
    # mean_type ∈ {:fixed_mean, :mean}
    # uncertainty_type ∈ {:fixed_variance, :variance, :precision, :log_variance}
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

function GaussianNode(; id=generateNodeId(GaussianNode), form::Symbol=:variance, m::Union{Float64,Vector{Float64}}=[NaN], V::Union{Float64,Matrix{Float64}}=reshape([NaN],1,1))
    mean_type = isValid(m) ? :fixed_mean : :mean
    if form == :variance
        uncertainty_type = isValid(V) ? :fixed_variance : :variance
    else
        uncertainty_type = form
    end

    num_interfaces = ((mean_type != :fixed_mean) ? 1 : 0) + ((uncertainty_type != :fixed_variance) ? 1 : 0) + 1
    self = GaussianNode{Val{mean_type},Val{uncertainty_type}}(id, Vector{Interface}(num_interfaces), Dict{Symbol,Interface}())
    addNode!(currentGraph(), self)

    next_interface_index = 1 # Counter keeping track of constructed interfaces

    if mean_type == :fixed_mean
        self.m = (typeof(m)==Float64) ? [m] : deepcopy(m)
    else
        self.i[:mean] = self.interfaces[next_interface_index] = Interface(self)
        next_interface_index += 1
    end

    if uncertainty_type == :fixed_variance
        self.V = (typeof(V)==Float64) ? reshape([V],1,1) : deepcopy(V)
    else
        self.i[uncertainty_type] = self.interfaces[next_interface_index] = Interface(self)
        next_interface_index += 1
    end

    self.i[:out] = self.interfaces[next_interface_index] = Interface(self)

    return self
end

isDeterministic(::GaussianNode) = false


############################################
# Sumproduct rules
############################################

"""
GaussianNode{:mean, :variance/:precision}:

     N       δ
    --->[N]<---
    <--  |
         v δ

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Any,
                            msg_var_prec::Message{DeltaDistribution{Float64}},
                            msg_out::Message{DeltaDistribution{Float64}})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var_prec.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Any,
                            msg_var_prec::Message{DeltaDistribution{Float64}},
                            msg_out::Message{DeltaDistribution{Float64}})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_var_prec.payload.m

    return outbound_dist
end

"""
GaussianNode{:mean, :variance}:

     δ       Ig
    --->[N]<---
         |  -->
         v δ

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::InverseGammaDistribution,
                            msg_mean::Message{DeltaDistribution{Float64}},
                            msg_var::Any,
                            msg_out::Message{DeltaDistribution{Float64}})

    y = msg_out.payload.m
    m = msg_mean.payload.m
    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

     δ      Gam
    --->[N]<---
         |  -->
         v δ

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GammaDistribution,
                            msg_mean::Message{DeltaDistribution{Float64}},
                            msg_prec::Any,
                            msg_out::Message{DeltaDistribution{Float64}})

    y = msg_out.payload.m
    m = msg_mean.payload.m
    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

"""
GaussianNode{:mean, :variance/:precision}:

     δ       δ
    --->[N]<---
         | |
       N v v

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Message{DeltaDistribution{Float64}},
                            msg_var_prec::Message{DeltaDistribution{Float64}},
                            msg_out::Any)

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var_prec.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Message{DeltaDistribution{Float64}},
                            msg_var_prec::Message{DeltaDistribution{Float64}},
                            msg_out::Any)

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_var_prec.payload.m

    return outbound_dist
end


#############################################
# Sum-product rules for fixed mean
#############################################

"""
GaussianNode{:fixed_mean, fixed_variance}:

        N
    [N]--->
        -->
"""
function sumProductRule!(   node::GaussianNode{Val{:fixed_mean},Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_out::Any)

    outbound_dist.m = node.m[1]
    outbound_dist.V = node.V[1,1]
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :variance}:

     Ig     δ
    --->[N]--->
    <--

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:fixed_mean},Val{:variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::InverseGammaDistribution,
                            msg_var::Any,
                            msg_out::Message{DeltaDistribution{Float64}})

    y = msg_out.payload.m
    m = node.m[1]
    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :precision}

    Gam     δ
    --->[N]--->
    <--

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:fixed_mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GammaDistribution,
                            msg_prec::Any,
                            msg_out::Message{DeltaDistribution{Float64}})

    y = msg_out.payload.m
    m = node.m[1]
    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(y-m)^2

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :variance/:precision}

     δ      N
    --->[N]--->
            -->

    Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; table 5.2
"""
function sumProductRule!(   node::GaussianNode{Val{:fixed_mean},Val{:variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_var::Message{DeltaDistribution{Float64}},
                            msg_out::Any)

    outbound_dist.m = node.m[1]
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!(   node::GaussianNode{Val{:fixed_mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_prec::Message{DeltaDistribution{Float64}},
                            msg_out::Any)

    outbound_dist.m = node.m[1]
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_prec.payload.m

    return outbound_dist
end


#############################################
# Sum-product rules for fixed variance
#############################################

"""
GaussianNode{:mean, :fixed_variance}

     N      δ
    --->[N]--->
    <--
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Any,
                            msg_out::Message{DeltaDistribution{Float64}})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussianDistribution{dims},
                                msg_mean::Any,
                                msg_out::Message{MvDeltaDistribution{Float64,dims}})

    outbound_dist.m[:] = msg_out.payload.m
    outbound_dist.V[:] = node.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

    return outbound_dist
end

"""
GaussianNode{:mean, :fixed_variance}

     N      N
    --->[N]--->
    <--
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Any,
                            msg_out::Message{GaussianDistribution})

    ensureParameters!(msg_out.payload, (:m,:V))
    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_out.payload.V + node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussianDistribution{dims},
                                msg_mean::Any,
                                msg_out::Message{MvGaussianDistribution{dims}})

    ensureParameters!(msg_out.payload, (:m,:V))
    outbound_dist.m[:] = msg_out.payload.m
    outbound_dist.V[:] = msg_out.payload.V + node.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

    return outbound_dist
end

"""
GaussianNode{:mean, :fixed_variance}

     δ      N
    --->[N]--->
            -->
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Message{DeltaDistribution{Float64}},
                            msg_out::Any)

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{2}},
                                outbound_dist::MvGaussianDistribution{dims},
                                msg_mean::Message{MvDeltaDistribution{Float64,dims}},
                                msg_out::Any)

    outbound_dist.m[:] = msg_mean.payload.m
    outbound_dist.V[:] = node.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

    return outbound_dist
end

"""
GaussianNode{:mean, :fixed_variance}

     N      N
    --->[N]--->
            -->
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            msg_mean::Message{GaussianDistribution},
                            msg_out::Any)

    ensureParameters!(msg_mean.payload, (:m,:V))
    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_mean.payload.V + node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{2}},
                                outbound_dist::MvGaussianDistribution{dims},
                                msg_mean::Message{MvGaussianDistribution{dims}},
                                msg_out::Any)

    ensureParameters!(msg_mean.payload, (:m,:V))
    outbound_dist.m[:] = msg_mean.payload.m
    outbound_dist.V[:] = msg_mean.payload.V + node.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

    return outbound_dist
end


############################################
# Naive variational update functions
############################################

"""
GaussianNode{:mean, :variance}:

     N     q(Ig)
    --->[N]<---
    <--  |
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::Any,
                            marg_variance::InverseGammaDistribution,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m,))

    outbound_dist.m = marg_out.m
    outbound_dist.V = (marg_variance.a + 1)/marg_variance.b
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

     N     q(Gam)
    --->[N]<---
    <--  |
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::Any,
                            marg_prec::GammaDistribution,
                            marg_y::GaussianDistribution)

    ensureParameters!(marg_y, (:m,))

    outbound_dist.m = marg_y.m
    outbound_dist.W = marg_prec.a / marg_prec.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

     N     q(W)
    --->[N]<---
    <--  |
         v q(N)
"""
function variationalRule!{dims}(node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussianDistribution{dims},
                                marg_mean::Any,
                                marg_prec::WishartDistribution{dims},
                                marg_y::MvGaussianDistribution{dims})

    ensureParameters!(marg_y, (:m,:W))
    outbound_dist.m = deepcopy(marg_y.m)
    outbound_dist.W = marg_prec.nu*marg_prec.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.V)

    return outbound_dist
end

"""
GaussianNode{:mean, :log_variance}:

     N     q(N)
    --->[N]<---
    <--  |
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:log_variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::Any,
                            marg_log_var::GaussianDistribution,
                            marg_y::GaussianDistribution)

    ensureParameters!(marg_y, (:m,))
    ensureParameters!(marg_log_var, (:m, :V))

    outbound_dist.m = marg_y.m
    V = exp(marg_log_var.m + 0.5*marg_log_var.V)
    outbound_dist.V = (V > huge ? huge : V)
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :variance}:

    q(N)    Ig
    --->[N]<---
         |  -->
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::InverseGammaDistribution,
                            marg_mean::GaussianDistribution,
                            marg_var::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :V))
    ensureParameters!(marg_mean, (:m, :V))

    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(marg_out.m - marg_mean.m)^2 + 0.5*(marg_mean.V + marg_out.V)

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

    q(N)    Gam
    --->[N]<---
         |  -->
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GammaDistribution,
                            marg_mean::GaussianDistribution,
                            marg_prec::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :V))
    ensureParameters!(marg_mean, (:m, :V))

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(marg_out.m - marg_mean.m)^2 + 0.5*(marg_mean.V + marg_out.V)

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

    q(N)     W
    --->[N]<---
         |  -->
         v q(N)
"""
function variationalRule!{dims}(node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{2}},
                                outbound_dist::WishartDistribution{dims},
                                marg_mean::MvGaussianDistribution{dims},
                                marg_prec::Any,
                                marg_out::MvGaussianDistribution{dims})

    ensureParameters!(marg_out, (:m, :V))
    ensureParameters!(marg_mean, (:m, :V))
    outbound_dist.nu = 2.0 + dimensions(marg_mean)
    outbound_dist.V = inv(cholfact( marg_out.V + marg_mean.V + (marg_out.m - marg_mean.m)*(marg_out.m - marg_mean.m)' ))

    return outbound_dist
end

"""
GaussianNode{:mean, :log_variance}:

    q(N)     N
    --->[N]<---
         |  -->
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:log_variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::GaussianDistribution,
                            marg_log_var::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :V))
    ensureParameters!(marg_mean, (:m, :V))

    outbound_dist.m = log(marg_mean.m^2 + marg_mean.V - 2.0*marg_mean.m*marg_out.m + marg_out.m^2 + marg_out.V)
    outbound_dist.V = 2.0
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :variance}:

    q(N)   q(Ig)
    --->[N]<---
         | |
       N v v
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::GaussianDistribution,
                            marg_var::InverseGammaDistribution,
                            marg_out::Any)

    ensureParameters!(marg_mean, (:m,))

    outbound_dist.m = marg_mean.m
    outbound_dist.V = marg_var.b / (marg_var.a-1.0)
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

    q(N)   q(G)
    --->[N]<---
         | |
       N v v
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::GaussianDistribution,
                            marg_prec::GammaDistribution,
                            marg_out::Any)

    ensureParameters!(marg_mean, (:m,))

    outbound_dist.m = marg_mean.m
    outbound_dist.W = marg_prec.a / marg_prec.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

    q(N)   q(W)
    --->[N]<---
         | |
       N v v
"""
function variationalRule!{dims}(node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{3}},
                                outbound_dist::MvGaussianDistribution{dims},
                                marg_mean::MvGaussianDistribution{dims},
                                marg_prec::WishartDistribution{dims},
                                marg_out::Any)

    ensureParameters!(marg_mean, (:m,:W))

    outbound_dist.m = deepcopy(marg_mean.m)
    outbound_dist.W = marg_prec.nu*marg_prec.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.V)

    return outbound_dist
end

"""
GaussianNode{:mean, :log_variance}:

    q(N)   q(N)
    --->[N]<---
         | |
       N v v
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:log_variance}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::GaussianDistribution,
                            marg_log_var::GaussianDistribution,
                            marg_out::Any)

    ensureParameters!(marg_mean, (:m,))
    ensureParameters!(marg_log_var, (:m,:V))

    outbound_dist.m = marg_mean.m
    V = exp(marg_log_var.m + 0.5*marg_log_var.V)
    outbound_dist.V = (V > huge ? huge : V)
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end


##########################################################
# Naive variational update functions with fixed interfaces
##########################################################

"""
GaussianNode{:mean, :fixed_variance}

     N     q(N)
    --->[N]--->
    <--
"""
function variationalRule!(  node::GaussianNode{Val{:mean}, Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m,))

    outbound_dist.m = marg_out.m
    outbound_dist.V = node.V[1,1]
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :variance}

     Ig    q(N)
    --->[N]--->
    <--
"""
function variationalRule!(  node::GaussianNode{Val{:fixed_mean},Val{:variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::InverseGammaDistribution,
                            marg_var::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :V))

    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(marg_out.m - node.m[1])^2 + 0.5*marg_out.V

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :precision}

    Gam    q(N)
    --->[N]--->
    <--
"""
function variationalRule!(  node::GaussianNode{Val{:fixed_mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GammaDistribution,
                            marg_prec::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :V))

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*(marg_out.m - node.m[1])^2 + 0.5*marg_out.V

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :log_variance}

     N     q(N)
    --->[N]--->
    <--
"""
function variationalRule!(  node::GaussianNode{Val{:fixed_mean},Val{:log_variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::GaussianDistribution,
                            marg_log_var::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :V))

    outbound_dist.m = log(node.m[1]^2 - 2.0*node.m[1]*marg_out.m + marg_out.m^2 + marg_out.V)
    outbound_dist.V = 2.0
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :fixed_variance}

    q(N)    N
    --->[N]--->
            -->
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            marg_mean::GaussianDistribution,
                            marg_out::Any)

    ensureParameters!(marg_mean, (:m,))

    outbound_dist.m = marg_mean.m
    outbound_dist.V = node.V[1,1]
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :variance}

    q(Ig)    N
     --->[N]--->
             -->
"""
function variationalRule!(  node::GaussianNode{Val{:fixed_mean},Val{:variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            marg_var::InverseGammaDistribution,
                            marg_out::Any)

    outbound_dist.m = node.m[1]
    outbound_dist.V = marg_var.b/(marg_var.a - 1.0)
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :precision}

    q(Gam)   N
     --->[N]--->
             -->
"""
function variationalRule!(  node::GaussianNode{Val{:fixed_mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            marg_prec::GammaDistribution,
                            marg_out::Any)

    outbound_dist.m = node.m[1]
    outbound_dist.V = NaN
    outbound_dist.xi = NaN
    outbound_dist.W = marg_prec.a / marg_prec.b

    return outbound_dist
end

"""
GaussianNode{:fixed_mean, :log_variance}

    q(N)    N
    --->[N]--->
            -->
"""
function variationalRule!(  node::GaussianNode{Val{:fixed_mean},Val{:log_variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            marg_log_var::GaussianDistribution,
                            marg_out::Any)

    ensureParameters!(marg_log_var, (:m,))

    outbound_dist.m = node.m[1]
    V = exp(marg_log_var.m + 0.5*marg_log_var.V)
    outbound_dist.V = (V > huge ? huge : V)
    outbound_dist.xi = NaN
    outbound_dist.W = NaN

    return outbound_dist
end


############################################
# Structured variational update functions
############################################

"""
GaussianNode{:mean, :precision}:

     St     Gam
    --->[N]<---
    <--  .
         .
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::StudentsTDistribution,
                            msg_mean::Any,
                            msg_prec::Message{GammaDistribution},
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :W))

    outbound_dist.m = marg_out.m
    outbound_dist.lambda = (2.0*marg_out.W*msg_prec.payload.a)/(2.0*marg_out.W*msg_prec.payload.b + 1)
    outbound_dist.nu = 2.0*msg_prec.payload.a

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

     N      Gam
    --->[N]<---
         .  -->
         .
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::GammaDistribution,
                            msg_mean::Message{GaussianDistribution},
                            msg_prec::Any,
                            marg_out::GaussianDistribution)

    ensureParameters!(marg_out, (:m, :W))
    ensureParameters!(msg_mean.payload, (:m,))

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*((1.0/marg_out.W) + (msg_mean.payload.m - marg_out.m)^2)

    return outbound_dist
end

"""
GaussianNode{:mean, precision}:

      q(NGam)
    --->[N]<---
         .
       N . |
         v v
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::GaussianDistribution,
                            marg::NormalGammaDistribution,
                            ::NormalGammaDistribution, # Same distribution as marg
                            marg_out::Any)

    outbound_dist.m = marg.m
    outbound_dist.W = marg.a / marg.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end