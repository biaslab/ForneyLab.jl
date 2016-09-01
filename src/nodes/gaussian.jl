export GaussianNode

"""
Description:

    Gaussian factor node.

Interfaces:

    The GaussianNode can take multiple parametrization forms, specified
    through the form argument:

    - form == :variance     => f(out,m,V)  = N(out|m,V)        (default)
    - form == :precision    => f(out,m,W)  = N(out|m,inv(W))

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
       - uncertainty_type ∈ {:fixed_variance, :variance, :precision}

Construction:

    GaussianNode(form=:variance, id=:my_node)           # N(out|m,V)        (default)
    GaussianNode(form=:precision, m=0.0, id=:my_node)   # N(out|0.0,inv(W)) (fixed mean)
    GaussianNode(V=1.0, id=:my_node)                    # N(out|m,1.0)      (fixed variance)
    GaussianNode(m=0.0, V=1.0, id=:my_node)             # N(out|0.0,1.0)    (fixed distribution)
"""
type GaussianNode{mean_type,uncertainty_type} <: Node
    # mean_type ∈ {:fixed_mean, :mean}
    # uncertainty_type ∈ {:fixed_variance, :variance, :precision}
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}
    # Fixed parameters
    m::Vector{Float64}
    V::AbstractMatrix{Float64}

    function GaussianNode(id::Symbol, interfaces::Array{Interface,1}, i::Dict{Symbol,Interface})
        new(id, interfaces, i)
    end
end

function GaussianNode(; id=generateNodeId(GaussianNode), form::Symbol=:variance, m::Union{Float64,Vector{Float64}}=[NaN], V::Union{Float64, AbstractMatrix{Float64}}=reshape([NaN],1,1))
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
                            outbound_dist::Gaussian,
                            msg_mean::Any,
                            msg_var::Message{Delta{Float64}},
                            msg_out::Message{Delta{Float64}})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            msg_mean::Any,
                            msg_prec::Message{Delta{Float64}},
                            msg_out::Message{Delta{Float64}})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_prec.payload.m

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Any,
                                msg_prec::Message{MatrixDelta{Float64,dims,dims}},
                                msg_out::Message{MvDelta{Float64,dims}})

    outbound_dist.m = deepcopy(msg_out.payload.m)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.V)
    outbound_dist.W = deepcopy(msg_prec.payload.M)

    return outbound_dist
end

"""
GaussianNode{:mean, :variance}

     N       δ
    --->[N]<---
    <--  | N
         v
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            msg_mean::Any,
                            msg_var::Message{Delta{Float64}},
                            msg_out::Message{Gaussian})

    ensureParameters!(msg_out.payload, (:m,:V))
    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_out.payload.V + msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}

     N       δ
    --->[N]<---
    <--  | N
         v
"""
function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Any,
                                msg_prec::Message{MatrixDelta{Float64, dims, dims}},
                                msg_out::Message{MvGaussian{dims}})

    ensureParameters!(msg_out.payload, (:m,:V))
    outbound_dist.m = deepcopy(msg_out.payload.m)
    outbound_dist.V = msg_out.payload.V + cholinv(msg_prec.payload.M)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

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
                            outbound_dist::InverseGamma,
                            msg_mean::Message{Delta{Float64}},
                            msg_var::Any,
                            msg_out::Message{Delta{Float64}})

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
                            outbound_dist::Gamma,
                            msg_mean::Message{Delta{Float64}},
                            msg_prec::Any,
                            msg_out::Message{Delta{Float64}})

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
                            outbound_dist::Gaussian,
                            msg_mean::Message{Delta{Float64}},
                            msg_var::Message{Delta{Float64}},
                            msg_out::Any)

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            msg_mean::Message{Delta{Float64}},
                            msg_prec::Message{Delta{Float64}},
                            msg_out::Any)

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = msg_prec.payload.m

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:variance}},
                                outbound_interface_index::Type{Val{3}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Message{MvDelta{Float64, dims}},
                                msg_var::Message{MatrixDelta{Float64, dims, dims}},
                                msg_out::Any)

    outbound_dist.m = deepcopy(msg_mean.payload.m)
    invalidate!(outbound_dist.xi)
    outbound_dist.V = deepcopy(msg_var.payload.M)
    invalidate!(outbound_dist.W)

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{3}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Message{MvDelta{Float64, dims}},
                                msg_prec::Message{MatrixDelta{Float64, dims, dims}},
                                msg_out::Any)

    outbound_dist.m = deepcopy(msg_mean.payload.m)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.V)
    outbound_dist.W = deepcopy(msg_prec.payload.M)

    return outbound_dist
end

"""
GaussianNode{:mean, :variance}

     N       δ
    --->[N]<---
       . | N
       v v
"""
function sumProductRule!(   node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            msg_mean::Message{Gaussian},
                            msg_var::Message{Delta{Float64}},
                            msg_out::Any)

    ensureParameters!(msg_mean.payload, (:m,:V))
    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_mean.payload.V + msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}

     N       δ
    --->[N]<---
       . | N
       v v
"""
function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{3}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Message{MvGaussian{dims}},
                                msg_prec::Message{MatrixDelta{Float64, dims, dims}},
                                msg_out::Any)

    ensureParameters!(msg_mean.payload, (:m,:V))
    outbound_dist.m = deepcopy(msg_mean.payload.m)
    outbound_dist.V = msg_mean.payload.V + cholinv(msg_prec.payload.M)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

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
                            outbound_dist::Gaussian,
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
                            outbound_dist::InverseGamma,
                            msg_var::Any,
                            msg_out::Message{Delta{Float64}})

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
                            outbound_dist::Gamma,
                            msg_prec::Any,
                            msg_out::Message{Delta{Float64}})

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
                            outbound_dist::Gaussian,
                            msg_var::Message{Delta{Float64}},
                            msg_out::Any)

    outbound_dist.m = node.m[1]
    outbound_dist.xi = NaN
    outbound_dist.V = msg_var.payload.m
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!(   node::GaussianNode{Val{:fixed_mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Gaussian,
                            msg_prec::Message{Delta{Float64}},
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
                            outbound_dist::Gaussian,
                            msg_mean::Any,
                            msg_out::Message{Delta{Float64}})

    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Any,
                                msg_out::Message{MvDelta{Float64,dims}})

    outbound_dist.m = deepcopy(msg_out.payload.m)
    outbound_dist.V = deepcopy(node.V)
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
                            outbound_dist::Gaussian,
                            msg_mean::Any,
                            msg_out::Message{Gaussian})

    ensureParameters!(msg_out.payload, (:m,:V))
    outbound_dist.m = msg_out.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = msg_out.payload.V + node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Any,
                                msg_out::Message{MvGaussian{dims}})

    ensureParameters!(msg_out.payload, (:m,:V))
    outbound_dist.m = deepcopy(msg_out.payload.m)
    outbound_dist.V = msg_out.payload.V + node.V
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
                            outbound_dist::Gaussian,
                            msg_mean::Message{Delta{Float64}},
                            msg_out::Any)

    outbound_dist.m = msg_mean.payload.m
    outbound_dist.xi = NaN
    outbound_dist.V = node.V[1,1]
    outbound_dist.W = NaN

    return outbound_dist
end

function sumProductRule!{dims}( node::GaussianNode{Val{:mean},Val{:fixed_variance}},
                                outbound_interface_index::Type{Val{2}},
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Message{MvDelta{Float64,dims}},
                                msg_out::Any)

    outbound_dist.m = deepcopy(msg_mean.payload.m)
    outbound_dist.V = deepcopy(node.V)
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
                            outbound_dist::Gaussian,
                            msg_mean::Message{Gaussian},
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
                                outbound_dist::MvGaussian{dims},
                                msg_mean::Message{MvGaussian{dims}},
                                msg_out::Any)

    ensureParameters!(msg_mean.payload, (:m,:V))
    outbound_dist.m = deepcopy(msg_mean.payload.m)
    outbound_dist.V = msg_mean.payload.V + node.V
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

    return outbound_dist
end


############################################
# Naive variational update functions
############################################

"""
GaussianNode{:mean, :precision}:

    q(N)   q(gam)
    --->[N]<---
         |
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            marg_mean::Any,
                            marg_prec::Univariate,
                            marg_out::Univariate)

    outbound_dist.m = unsafeMean(marg_out)
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = unsafeMean(marg_prec)

    return outbound_dist
end

function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Gamma,
                            marg_mean::Univariate,
                            marg_prec::Any,
                            marg_out::Univariate)

    outbound_dist.a = 1.5
    outbound_dist.b = 0.5*( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2 )

    return outbound_dist
end

function variationalRule!(  node::GaussianNode{Val{:mean},Val{:precision}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            marg_mean::Univariate,
                            marg_prec::Univariate,
                            marg_out::Any)

    outbound_dist.m = unsafeMean(marg_mean)
    outbound_dist.xi = NaN
    outbound_dist.V = NaN
    outbound_dist.W = unsafeMean(marg_prec)

    return outbound_dist
end

"""
GaussianNode{:mean, :variance}:

    q(N)   q(Ig)
    --->[N]<---
         |
         v q(N)
"""
function variationalRule!(  node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Gaussian,
                            marg_mean::Any,
                            marg_var::Univariate,
                            marg_out::Univariate)

    outbound_dist.m = unsafeMean(marg_out)
    outbound_dist.xi = NaN
    outbound_dist.V = unsafeMean(marg_var)
    outbound_dist.W = NaN

    return outbound_dist
end

function variationalRule!(  node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::InverseGamma,
                            marg_mean::Univariate,
                            marg_var::Any,
                            marg_out::Univariate)

    outbound_dist.a = -0.5
    outbound_dist.b = 0.5*(unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))^2)

    return outbound_dist
end

function variationalRule!(  node::GaussianNode{Val{:mean},Val{:variance}},
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            marg_mean::Univariate,
                            marg_var::Univariate,
                            marg_out::Any)

    outbound_dist.m = unsafeMean(marg_mean)
    outbound_dist.xi = NaN
    outbound_dist.V = unsafeMean(marg_var)
    outbound_dist.W = NaN

    return outbound_dist
end

"""
GaussianNode{:mean, :precision}:

    q(N)   q(W)
    --->[N]<---
         |
         v q(N)
"""
function variationalRule!{dims}(node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{1}},
                                outbound_dist::MvGaussian{dims},
                                marg_mean::Any,
                                marg_prec::MatrixVariate{dims, dims},
                                marg_out::Multivariate{dims})

    outbound_dist.m = unsafeMean(marg_out)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.V)
    outbound_dist.W = unsafeMean(marg_prec)

    return outbound_dist
end

function variationalRule!{dims}(node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{2}},
                                outbound_dist::Wishart{dims},
                                marg_mean::Multivariate{dims},
                                marg_prec::Any,
                                marg_out::Multivariate{dims})

    outbound_dist.nu = dims + 2.0
    outbound_dist.V = cholinv( unsafeCov(marg_out) + unsafeCov(marg_mean) + (unsafeMean(marg_out) - unsafeMean(marg_mean))*(unsafeMean(marg_out) - unsafeMean(marg_mean))' )

    return outbound_dist
end

function variationalRule!{dims}(node::GaussianNode{Val{:mean},Val{:precision}},
                                outbound_interface_index::Type{Val{3}},
                                outbound_dist::MvGaussian{dims},
                                marg_mean::Multivariate{dims},
                                marg_prec::MatrixVariate{dims, dims},
                                marg_out::Any)

    outbound_dist.m = unsafeMean(marg_mean)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.V)
    outbound_dist.W = unsafeMean(marg_prec)

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
                            outbound_dist::StudentsT,
                            msg_mean::Any,
                            msg_prec::Message{Gamma},
                            marg_out::Gaussian)

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
                            outbound_dist::Gamma,
                            msg_mean::Message{Gaussian},
                            msg_prec::Any,
                            marg_out::Gaussian)

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
                            outbound_dist::Gaussian,
                            marg::NormalGamma,
                            ::NormalGamma, # Same distribution as marg
                            marg_out::Any)

    outbound_dist.m = marg.m
    outbound_dist.W = marg.a / marg.b
    outbound_dist.xi = NaN
    outbound_dist.V = NaN

    return outbound_dist
end


############################
# Average energy functionals
############################

function U(::GaussianNode{Val{:mean},Val{:precision}}, dist_mean::Gaussian, dist_prec::Gamma, dist_y::Gaussian)
    ensureParameters!(dist_mean, (:m, :V))
    ensureParameters!(dist_y, (:m, :V))

    return  (1/2)*log(2*pi) -
            0.5*digamma(dist_prec.a) +
            0.5*log(dist_prec.b) +
            0.5*(dist_prec.a/dist_prec.b)*( dist_y.V + dist_mean.V + (dist_y.m - dist_mean.m)^2)
end

function U{dims}(::GaussianNode{Val{:mean},Val{:precision}}, dist_mean::MvGaussian{dims}, dist_prec::Wishart{dims}, dist_y::MvGaussian{dims})
    ensureParameters!(dist_mean, (:m, :V))
    ensureParameters!(dist_y, (:m, :V))

    return  (dims/2)*log(2*pi) -
            0.5*sum([digamma(0.5*(dist_prec.nu + 1 - i)) for i = 1:dims]) -
            (dims/2)*log(2) -
            0.5*log(det(dist_prec.V)) +
            0.5*trace( dist_prec.nu*dist_prec.V*( dist_y.V + dist_mean.V + (dist_y.m - dist_mean.m)*(dist_y.m - dist_mean.m)' ) )
end