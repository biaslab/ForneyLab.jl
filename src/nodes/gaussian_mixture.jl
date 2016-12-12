export GaussianMixtureNode

"""
Description:
    Gaussian mixture model
    Node function for two clusters:
    f(m,w,pi,x,z) = N(x|m1,1/w1)^z[1] * N(x|m2,1/w2)^z[2]


          _________________
      m   |               |  x
     -----|               |----
          |               |  z
          |       GM      |----
          |               |
          |               |
      w   |               |
     -----|               |
          |_______________|
Interfaces:
    1 i[:m]
    2 i[:w]
    3 i[:x]
    4 i[:z]
Construction:
    GaussianMixtureNode(id=:my_node)
"""
type GaussianMixtureNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function GaussianMixtureNode(; id=generateNodeId(GaussianMixtureNode))
        self = new(id, Array(Interface, 4), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:m, :w, :x, :z])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

ForneyLab.isDeterministic(::GaussianMixtureNode) = false

############################################
# Sumproduct rules
############################################
#Sum product rule
#Predictive distribution towards x
#Univariate case with 2 clusters
function sumProductRule!(   node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Mixture{Gaussian},
                            msg_m::Message{Partitioned{Gaussian,2}},
                            msg_w::Message{Partitioned{Gamma,2}},
                            msg_x::Any,
                            msg_z::Message{Bernoulli})

    #Ensure that the messages are the right form
    ensureParameters!(msg_m.payload.factors[1], (:m, :V))
    ensureParameters!(msg_m.payload.factors[2], (:m, :V))
    resize!(outbound_dist, 2)

    #Calculate the mean, variance and weight for the first component
    outbound_dist.components[1].m   = unsafeMean(msg_m.payload.factors[1])
    outbound_dist.components[1].V   = msg_w.payload.factors[1].b/(msg_w.payload.factors[1].a-1) + unsafeCov(msg_m.payload.factors[1])
    outbound_dist.components[1].xi  = NaN
    outbound_dist.components[1].W   = NaN

    outbound_dist.weights[1]      = unsafeMean(msg_z.payload)

    #Calculate the mean, variance and weight for the second component
    outbound_dist.components[2].m   = unsafeMean(msg_m.payload.factors[2])
    outbound_dist.components[2].V   = msg_w.payload.factors[2].b/(msg_w.payload.factors[2].a - 1) + unsafeCov(msg_m.payload.factors[2])
    outbound_dist.components[2].xi  = NaN
    outbound_dist.components[2].W   = NaN

    outbound_dist.weights[2]      = 1 - unsafeMean(msg_z.payload)

    return outbound_dist
end

#Sum product predictive distribution
#Univariate case with k clusters
function sumProductRule!{n_factors}(   node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Mixture{Gaussian},
                            msg_m::Message{Partitioned{Gaussian,n_factors}},
                            msg_w::Message{Partitioned{Gamma,n_factors}},
                            msg_x::Any,
                            msg_z::Message{Categorical{n_factors}})

    #Ensure that the outbound distribution has the correct form
    resize!(outbound_dist, n_factors)
    w = ones(n_factors)

    #Calculate the mean and variance for each component
    for k = 1:n_factors
        outbound_dist.components[k].m   = unsafeMean(msg_m.payload.factors[k])
        outbound_dist.components[k].V   = msg_w.payload.factors[k].b/(msg_w.payload.factors[k].a - 1) + unsafeCov(msg_m.payload.factors[k])
        outbound_dist.components[k].xi  = NaN
        outbound_dist.components[k].W   = NaN

        outbound_dist.weights[k]        = unsafeMean(msg_z.payload)[k]

    end

    return outbound_dist
end

#Sum product predictive distribution
#Multivariate case with 2 clusters
function sumProductRule!{n_factors,dims}(   node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Mixture{MvGaussian{dims}},
                            msg_m::Message{Partitioned{MvGaussian{dims},n_factors}},
                            msg_w::Message{Partitioned{Wishart{dims},n_factors}},
                            msg_x::Any,
                            msg_z::Message{Bernoulli})

    #Ensure that the parameters are in the right form
    resize!(outbound_dist, 2)
    w  =  ones(2)

    #Calculate the mean and variance for the first component
    outbound_dist.components[1].m = unsafeMean(msg_m.payload.factors[1])
    outbound_dist.components[1].V = inv(msg_w.payload.factors[1].V)/(msg_w.payload.factors[1].nu - dims - 1.) + unsafeCov(msg_m.payload.factors[1])
    invalidate!(outbound_dist.components[1].xi)
    invalidate!(outbound_dist.components[1].W)

    #Calculate the mean and  variance for the second component
    outbound_dist.components[2].m = unsafeMean(msg_m.payload.factors[2])
    outbound_dist.components[2].V = inv(msg_w.payload.factors[2].V)/(msg_w.payload.factors[2].nu - dims - 1.) + unsafeCov(msg_m.payload.factors[2])
    invalidate!(outbound_dist.components[2].xi)
    invalidate!(outbound_dist.components[2].W)

    #Calculating the weights
    outbound_dist.weights[1]     = unsafeMean(msg_z.payload)
    outbound_dist.weights[2]     = 1 - unsafeMean(msg_z.payload)

    return outbound_dist
end

#Sum product predictive distribution
#Multivariate case with k clusters
function sumProductRule!{n_factors,dims}(   node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Mixture{MvGaussian{dims}},
                            msg_m::Message{Partitioned{MvGaussian{dims},n_factors}},
                            msg_w::Message{Partitioned{Wishart{dims},n_factors}},
                            msg_x::Any,
                            msg_z::Message{Categorical{n_factors}})

    #Ensure outbound message has the right form
    resize!(outbound_dist, n_factors)
    w = ones(n_factors)

    #Calculate the mean, variance and weights of all components
    for k = 1:n_factors
        outbound_dist.components[k].m = unsafeMean(msg_m.payload.factors[k])
        outbound_dist.components[k].V = inv(msg_w.payload.factors[k].V)/(msg_w.payload.factors[k].nu - dims - 1.) + unsafeCov(msg_m.payload.factors[k])
        invalidate!(outbound_dist.components[k].xi)
        invalidate!(outbound_dist.components[k].W)

        outbound_dist.weights[k]     = unsafeMean(msg_z.payload)[k]
    end

    return outbound_dist
end

# ############################################
# # Naive variational update functions
# ############################################

function GMBackwardMRule!(outbound_dist::Gaussian, q_w_k::Univariate, q_x::Univariate, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    outbound_dist.m  = unsafeMean(q_x)
    outbound_dist.W  = z_k_hat*unsafeMean(q_w_k)
    outbound_dist.V  = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end

function GMBackwardMRule!{dims}(outbound_dist::MvGaussian{dims}, q_w_k::MatrixVariate{dims, dims}, q_x::Multivariate{dims}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    outbound_dist.m = unsafeMean(q_x)
    outbound_dist.W = z_k_hat*unsafeMean(q_w_k)
    invalidate!(outbound_dist.V)
    invalidate!(outbound_dist.xi)

    return outbound_dist
end


# VMP message towards i[:m]
# Univariate gaussian with two clusters
function variationalRule!{n_factors, T<:Univariate}(node::GaussianMixtureNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::Partitioned{Gaussian, n_factors},
                                                    q_m::Any,
                                                    q_w::Partitioned{T, n_factors},
                                                    q_x::Univariate,
                                                    q_z::Bernoulli)

    GMBackwardMRule!(outbound_dist.factors[1], q_w.factors[1], q_x, unsafeMean(q_z))
    GMBackwardMRule!(outbound_dist.factors[2], q_w.factors[2], q_x, 1.0 - unsafeMean(q_z))

    return outbound_dist
end

# VMP message towards i[:m]
# Multivariate Gaussian with two clusters
# TODO: constrain to MatrixVariate{dims, dims}
function variationalRule!{dims, n_factors, T<:MatrixVariate}(   node::GaussianMixtureNode,
                                                                outbound_interface_index::Type{Val{1}},
                                                                outbound_dist::Partitioned{MvGaussian{dims}, n_factors},
                                                                q_m::Any,
                                                                q_w::Partitioned{T, n_factors},
                                                                q_x::Multivariate{dims},
                                                                q_z::Bernoulli)

    GMBackwardMRule!(outbound_dist.factors[1], q_w.factors[1], q_x, unsafeMean(q_z))
    GMBackwardMRule!(outbound_dist.factors[2], q_w.factors[2], q_x, 1.0 - unsafeMean(q_z))

    return outbound_dist
end
#
#
#
#
# VMP message towards i[:m]
# Multivariate Gaussian with multiple clusters
# TODO: constrain to MatrixVariate{dims, dims}
function variationalRule!{dims, n_factors, T<:MatrixVariate}(   node::GaussianMixtureNode,
                                                                outbound_interface_index::Type{Val{1}},
                                                                outbound_dist::Partitioned{MvGaussian{dims}, n_factors},
                                                                q_m::Any,
                                                                q_w::Partitioned{T, n_factors},
                                                                q_x::Multivariate{dims},
                                                                q_z::Categorical{n_factors})
    z_mean = unsafeMean(q_z)
    for k = 1:n_factors
        GMBackwardMRule!(outbound_dist.factors[k], q_w.factors[k], q_x, z_mean[k])
    end

    return outbound_dist
end

# VMP message towards i[:m]
# Univariate Gaussian with multiple clusters
function variationalRule!{n_factors, T<:Univariate}(node::GaussianMixtureNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::Partitioned{Gaussian, n_factors},
                                                    q_m::Any,
                                                    q_w::Partitioned{T, n_factors},
                                                    q_x::Univariate,
                                                    q_z::Categorical{n_factors})
    z_mean = unsafeMean(q_z)
    for k=1:n_factors
        GMBackwardMRule!(outbound_dist.factors[k], q_w.factors[k], q_x, z_mean[k])
    end

    return outbound_dist
end


function GMBackwardWRule!(outbound_dist::Gamma, q_m_k::Univariate, q_x::Univariate, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    outbound_dist.a   =   1.0 + 0.5*z_k_hat
    outbound_dist.b   =   0.5*z_k_hat*( (unsafeMean(q_x) - unsafeMean(q_m_k))^2 + unsafeCov(q_x) + unsafeCov(q_m_k) )

    return outbound_dist
end

function GMBackwardWRule!{dims}(outbound_dist::Wishart{dims}, q_m_k::Multivariate{dims}, q_x::Multivariate{dims}, z_k_hat::Float64)
    z_k_hat = clamp(z_k_hat, tiny, 1.0 - tiny)

    outbound_dist.nu = 1.0 + z_k_hat + dims
    outbound_dist.V  = cholinv( z_k_hat*( (unsafeMean(q_x) - unsafeMean(q_m_k))*(unsafeMean(q_x) - unsafeMean(q_m_k))' + unsafeCov(q_x) + unsafeCov(q_m_k) ) )

    return outbound_dist
end


# VMP message towards i[:w]
# Univariate gaussian with two clusters
function variationalRule!{n_factors, T<:Univariate}(node::GaussianMixtureNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::Partitioned{Gamma, n_factors},
                                                    q_m::Partitioned{T, n_factors},
                                                    q_w::Any,
                                                    q_x::Univariate,
                                                    q_z::Bernoulli)

    GMBackwardWRule!(outbound_dist.factors[1], q_m.factors[1], q_x, unsafeMean(q_z))
    GMBackwardWRule!(outbound_dist.factors[2], q_m.factors[2], q_x, 1.0 - unsafeMean(q_z))

    return outbound_dist
end
#
#
#
# VMP message towards i[:w]
# Multivariate Gaussian with two clusters
# TODO: constrain to Multivariate{dims}
function variationalRule!{dims, n_factors, T<:Multivariate}(node::GaussianMixtureNode,
                                                            outbound_interface_index::Type{Val{2}},
                                                            outbound_dist::Partitioned{Wishart{dims}, n_factors},
                                                            q_m::Partitioned{T, n_factors},
                                                            q_w::Any,
                                                            q_x::Multivariate{dims},
                                                            q_z::Bernoulli)

    GMBackwardWRule!(outbound_dist.factors[1], q_m.factors[1], q_x, unsafeMean(q_z))
    GMBackwardWRule!(outbound_dist.factors[2], q_m.factors[2], q_x, 1.0 - unsafeMean(q_z))

    return outbound_dist
end

#
#
# VMP message towards i[:w]
# Multivariate Gaussian with multiple clusters
# TODO: constrain to Multivariate{dims}
function variationalRule!{dims, n_factors, T<:Multivariate}(node::GaussianMixtureNode,
                                                            outbound_interface_index::Type{Val{2}},
                                                            outbound_dist::Partitioned{Wishart{dims}, n_factors},
                                                            q_m::Partitioned{T, n_factors},
                                                            q_w::Any,
                                                            q_x::Multivariate{dims},
                                                            q_z::Categorical{n_factors})
    z_mean = unsafeMean(q_z)
    for k = 1:n_factors
        GMBackwardWRule!(outbound_dist.factors[k], q_m.factors[k], q_x, z_mean[k])
    end

    return outbound_dist
end

# VMP message towards i[:w]
# Univariate Gaussian with multiple clusters
function variationalRule!{n_factors, T<:Univariate}(node::GaussianMixtureNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::Partitioned{Gamma, n_factors},
                                                    q_m::Partitioned{T, n_factors},
                                                    q_w::Any,
                                                    q_x::Univariate,
                                                    q_z::Categorical{n_factors})
    z_mean = unsafeMean(q_z)
    for k = 1:n_factors
        GMBackwardWRule!(outbound_dist.factors[k], q_m.factors[k], q_x, z_mean[k])
    end

    return outbound_dist
end


function GMBackwardZRule!(outbound_dist::Bernoulli, q_m::Partitioned, q_w::Partitioned, q_x::ProbabilityDistribution)
    rho = zeros(2)
    for k = 1:2
        rho[k] = clamp(exp(-averageEnergy(GaussianNode, q_m.factors[k], q_w.factors[k], q_x)), tiny, huge)
    end
    outbound_dist.p = rho[1]/sum(rho)

    return outbound_dist
end

function GMBackwardZRule!{n_factors}(outbound_dist::Categorical{n_factors}, q_m::Partitioned, q_w::Partitioned, q_x::ProbabilityDistribution)
    rho = zeros(n_factors)
    for k = 1:n_factors
        rho[k] = clamp(exp(-averageEnergy(GaussianNode, q_m.factors[k], q_w.factors[k], q_x)), tiny, huge)
    end
    outbound_dist.p = rho./sum(rho)

    return outbound_dist
end


# VMP message towards i[:z]
# Univariate gaussian with two clusters
function variationalRule!{n_factors, TN<:Univariate, TG<:Univariate}(   node::GaussianMixtureNode,
                                                                        outbound_interface_index::Type{Val{4}},
                                                                        outbound_dist::Bernoulli,
                                                                        q_m::Partitioned{TN, n_factors},
                                                                        q_w::Partitioned{TG, n_factors},
                                                                        q_x::Univariate,
                                                                        q_z::Any)
    GMBackwardZRule!(outbound_dist, q_m, q_w, q_x)

    return outbound_dist
end

# VMP message towards i[:z]
# Multivariate Gaussian with two clusters
# TODO: constrain to Multivariate{dims}, MatrixVariate{dims, dims}
function variationalRule!{dims, n_factors, TN<:Multivariate, TW<:MatrixVariate}(node::GaussianMixtureNode,
                                                                                outbound_interface_index::Type{Val{4}},
                                                                                outbound_dist::Bernoulli,
                                                                                q_m::Partitioned{TN,n_factors},
                                                                                q_w::Partitioned{TW, n_factors},
                                                                                q_x::Multivariate{dims},
                                                                                q_z::Any)
    GMBackwardZRule!(outbound_dist, q_m, q_w, q_x)

    return outbound_dist
end

# VMP message towards i[:z]
# Multivariate Gaussian with multiple clusters
# TODO: constrain to Multivariate{dims}, MatrixVariate{dims, dims}
function variationalRule!{dims, n_factors, TN<:Multivariate, TW<:MatrixVariate}(node::GaussianMixtureNode,
                                                                                outbound_interface_index::Type{Val{4}},
                                                                                outbound_dist::Categorical{n_factors},
                                                                                q_m::Partitioned{TN, n_factors},
                                                                                q_w::Partitioned{TW, n_factors},
                                                                                q_x::MvGaussian{dims},
                                                                                q_z::Any)
    GMBackwardZRule!(outbound_dist, q_m, q_w, q_x)

    return outbound_dist
end

# VMP message towards i[:z]
# Univariate gaussian with multiple clusters
function variationalRule!{n_factors, TN<:Univariate, TG<:Univariate}(   node::GaussianMixtureNode,
                                                                        outbound_interface_index::Type{Val{4}},
                                                                        outbound_dist::Categorical{n_factors},
                                                                        q_m::Partitioned{TN, n_factors},
                                                                        q_w::Partitioned{TG, n_factors},
                                                                        q_x::Univariate,
                                                                        q_z::Any)
    GMBackwardZRule!(outbound_dist, q_m, q_w, q_x)

    return outbound_dist
end


function GMForwardXRule!(outbound_dist::Gaussian, q_m::Partitioned, q_w::Partitioned, z_hat::Vector{Float64})
    n_factors = length(z_hat)
    W  = 0.0
    xi = 0.0
    for k = 1:n_factors
        W  += z_hat[k]*unsafeMean(q_w.factors[k])
        xi += unsafeMean(q_w.factors[k])*unsafeMean(q_m.factors[k])*z_hat[k]
    end

    outbound_dist.m  = NaN
    outbound_dist.V  = NaN
    outbound_dist.xi = xi
    outbound_dist.W  = W

    return outbound_dist
end

function GMForwardXRule!{dims}(outbound_dist::MvGaussian{dims}, q_m::Partitioned, q_w::Partitioned, z_hat::Vector{Float64})
    n_factors = length(z_hat)
    W  = Diagonal(zeros(dims))
    xi = zeros(dims)
    for k = 1:n_factors
        W  += z_hat[k]*unsafeMean(q_w.factors[k])
        xi += unsafeMean(q_w.factors[k])*unsafeMean(q_m.factors[k])*z_hat[k]
    end

    invalidate!(outbound_dist.m)
    invalidate!(outbound_dist.V)
    outbound_dist.xi = xi
    outbound_dist.W  = W

    return outbound_dist
end


# VMP message towards i[:x]
function variationalRule!{n_factors, TN<:Univariate, TG<:Univariate}(   node::GaussianMixtureNode,
                                                                        outbound_interface_index::Type{Val{3}},
                                                                        outbound_dist::Gaussian,
                                                                        q_m::Partitioned{TN, n_factors},
                                                                        q_w::Partitioned{TG, n_factors},
                                                                        q_x::Any,
                                                                        q_z::Bernoulli)

    GMForwardXRule!(outbound_dist, q_m, q_w, [unsafeMean(q_z), 1.0 - unsafeMean(q_z)])

    return outbound_dist
end
#
# #
# VMP message towards i[:x]
# TODO: replace MvGaussian{dims} with TN
function variationalRule!{dims, n_factors, TW<:MatrixVariate}(  node::GaussianMixtureNode,
                                                                outbound_interface_index::Type{Val{3}},
                                                                outbound_dist::MvGaussian{dims},
                                                                q_m::Partitioned{MvGaussian{dims}, n_factors},
                                                                q_w::Partitioned{TW, n_factors},
                                                                q_x::Any,
                                                                q_z::Bernoulli)

    GMForwardXRule!(outbound_dist, q_m, q_w, [unsafeMean(q_z), 1.0 - unsafeMean(q_z)])

    return outbound_dist
end
#
#
#
#
# #
# VMP message towards i[:x]
# TODO: replace MvGaussian{dims} with TN
function variationalRule!{dims, n_factors, TW<:MatrixVariate}(  node::GaussianMixtureNode,
                                                                outbound_interface_index::Type{Val{3}},
                                                                outbound_dist::MvGaussian{dims},
                                                                q_m::Partitioned{MvGaussian{dims}, n_factors},
                                                                q_w::Partitioned{TW, n_factors},
                                                                q_x::Any,
                                                                q_z::Categorical{n_factors})

    GMForwardXRule!(outbound_dist, q_m, q_w, unsafeMean(q_z))

    return outbound_dist
end
#
# VMP message towards i[:x]
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            q_m::Partitioned{Gaussian,n_factors},
                            q_w::Partitioned{Gamma,n_factors},
                            q_x::Any,
                            q_z::Categorical{n_factors})

    GMForwardXRule!(outbound_dist, q_m, q_w, unsafeMean(q_z))

    return outbound_dist
end
