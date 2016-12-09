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
# function sumProductRule!(   node::GaussianMixtureNode,
#                             outbound_interface_index::Type{Val{4}},
#                             outbound_dist::Mixture{Gaussian},
#                             msg_pi::Message{Beta},
#                             msg_m::Message{Partitioned{Gaussian,2}},
#                             msg_w::Message{Partitioned{Gamma,2}},
#                             msg_x::Any,
#                             msg_z::Message{Bernoulli})
#
#     #Ensure that the messages are the right form
#     ensureParameters!(msg_m.payload.factors[1], (:m, :V))
#     ensureParameters!(msg_m.payload.factors[2], (:m, :V))
#     resize!(outbound_dist, 2)
#
#     #Calculate the mean, variance and weight for the first component
#     outbound_dist.components[1].m   = msg_m.payload.factors[1].m
#     outbound_dist.components[1].V   = (msg_w.payload.factors[1].b/(msg_w.payload.factors[1].a-1) + msg_m.payload.factors[1].V)
#     outbound_dist.components[1].xi  = NaN
#     outbound_dist.components[1].W   = NaN
#
#     w1     = msg_z.payload.p*msg_pi.payload.a/(msg_pi.payload.a+msg_pi.payload.b)
#
#     #Calculate the mean, variance and weight for the second component
#     outbound_dist.components[2].m   = msg_m.payload.factors[2].m
#     outbound_dist.components[2].V   = (msg_w.payload.factors[2].b/(msg_w.payload.factors[2].a - 1) + msg_m.payload.factors[2].V)
#     outbound_dist.components[2].xi  = NaN
#     outbound_dist.components[2].W   = NaN
#
#     w2      = (1 - msg_z.payload.p)*msg_pi.payload.b/(msg_pi.payload.a + msg_pi.payload.b)
#
#
#     #Normalize the weights
#     outbound_dist.weights[1]        = w1/(w1 + w2)
#     outbound_dist.weights[2]        = w2/(w1 + w2)
#
#     return outbound_dist
# end
#
# #Sum product predictive distribution
# #Univariate case with k clusters
# function sumProductRule!{n_factors}(   node::GaussianMixtureNode,
#                             outbound_interface_index::Type{Val{4}},
#                             outbound_dist::Mixture{Gaussian},
#                             msg_pi::Message{Dirichlet{n_factors}},
#                             msg_m::Message{Partitioned{Gaussian,n_factors}},
#                             msg_w::Message{Partitioned{Gamma,n_factors}},
#                             msg_x::Any,
#                             msg_z::Message{Categorical{n_factors}})
#
#     #Ensure that the outbound distribution has the correct form
#     resize!(outbound_dist, n_factors)
#     w = ones(n_factors)
#
#     k = collect(1:n_factors)
#     sum_a = sum(msg_pi.payload.alpha[k])
#
#     #Calculate the mean and variance for each component
#     for k = 1:n_factors
#
#         ensureParameters!(msg_m.payload.factors[k], (:m, :V))
#         outbound_dist.components[k].m   = msg_m.payload.factors[k].m
#         outbound_dist.components[k].V   = (msg_w.payload.factors[k].b/(msg_w.payload.factors[k].a - 1) + msg_m.payload.factors[k].V)
#         outbound_dist.components[k].xi  = NaN
#         outbound_dist.components[k].W   = NaN
#
#         w[k]     = msg_z.payload.p[k]*msg_pi.payload.alpha[k]/(sum_a)
#
#     end
#
#     sum_w = sum(w)
#
#     #Normalize the weights
#     for k = 1:n_factors
#         outbound_dist.weights[k]        = w[k]/(sum_w)
#     end
#
#
#     return outbound_dist
# end
#
# #Sum product predictive distribution
# #Multivariate case with 2 clusters
# function sumProductRule!{n_factors,dims}(   node::GaussianMixtureNode,
#                             outbound_interface_index::Type{Val{4}},
#                             outbound_dist::Mixture{MvGaussian{dims}},
#                             msg_pi::Message{Beta},
#                             msg_m::Message{Partitioned{MvGaussian{dims},n_factors}},
#                             msg_w::Message{Partitioned{Wishart{dims},n_factors}},
#                             msg_x::Any,
#                             msg_z::Message{Bernoulli})
#
#     #Ensure that the parameters are in the right form
#     resize!(outbound_dist, 2)
#     ensureParameters!(msg_m.payload.factors[1], (:m, :V))
#     ensureParameters!(msg_m.payload.factors[2], (:m, :V))
#
#     w  =  ones(2)
#
#     #Calculate the mean and variance for the first component
#     outbound_dist.components[1].m = msg_m.payload.factors[1].m
#     outbound_dist.components[1].V = inv(msg_w.payload.factors[1].V)/(msg_w.payload.factors[1].nu - dims - 1.) + msg_m.payload.factors[1].V
#     invalidate!(outbound_dist.components[1].xi)
#     invalidate!(outbound_dist.components[1].W)
#
#     #Calculate the mean and  variance for the second component
#     outbound_dist.components[2].m = msg_m.payload.factors[2].m
#     outbound_dist.components[2].V = inv(msg_w.payload.factors[2].V)/(msg_w.payload.factors[2].nu - dims - 1.) + msg_m.payload.factors[2].V
#     invalidate!(outbound_dist.components[2].xi)
#     invalidate!(outbound_dist.components[2].W)
#
#     #Calculating the weights
#     w[1]     = msg_z.payload.p*msg_pi.payload.a/(msg_pi.payload.a + msg_pi.payload.b)
#     w[2]     = (1 - msg_z.payload.p)*msg_pi.payload.b/(msg_pi.payload.a + msg_pi.payload.b)
#
#
#     #Normalize the weights
#     outbound_dist.weights[1]  = w[1]/(w[1] + w[2])
#     outbound_dist.weights[2]  = w[2]/(w[1] + w[2])
#
#
#     return outbound_dist
# end
#
# #Sum product predictive distribution
# #Multivariate case with k clusters
# function sumProductRule!{n_factors,dims}(   node::GaussianMixtureNode,
#                             outbound_interface_index::Type{Val{4}},
#                             outbound_dist::Mixture{MvGaussian{dims}},
#                             msg_pi::Message{Dirichlet{n_factors}},
#                             msg_m::Message{Partitioned{MvGaussian{dims},n_factors}},
#                             msg_w::Message{Partitioned{Wishart{dims},n_factors}},
#                             msg_x::Any,
#                             msg_z::Message{Categorical{n_factors}})
#
#     #Ensure outbound message has the right form
#     resize!(outbound_dist, n_factors)
#     w = ones(n_factors)
#
#     k = collect(1:n_factors)
#     sum_a = sum(msg_pi.payload.alpha[k])
#
#     #Calculate the mean, variance and weights of all components
#     for k = 1:n_factors
#         ensureParameters!(msg_m.payload.factors[k], (:m, :V))
#         outbound_dist.components[k].m = msg_m.payload.factors[k].m
#         outbound_dist.components[k].V = inv(msg_w.payload.factors[k].V)/(msg_w.payload.factors[k].nu - dims - 1.) + msg_m.payload.factors[k].V
#         invalidate!(outbound_dist.components[k].xi)
#         invalidate!(outbound_dist.components[k].W)
#
#         w[k]     = msg_z.payload.p[k]*msg_pi.payload.alpha[k]/(sum_a)
#     end
#
#     sum_w = sum(w)
#
#     #Normalize the weights
#     for k = 1:n_factors
#         outbound_dist.weights[k]  = w[k]/(sum_w)
#     end
#
#
#     return outbound_dist
# end
#
# ############################################
# # Naive variational update functions
# ############################################
# VMP message towards i[:m]
# Univariate gaussian with two clusters
function variationalRule!{n_factors, T<:Univariate}(node::GaussianMixtureNode,
                                                    outbound_interface_index::Type{Val{1}},
                                                    outbound_dist::Partitioned{Gaussian, n_factors},
                                                    q_m::Any,
                                                    q_w::Partitioned{T, n_factors},
                                                    q_x::Univariate,
                                                    q_z::Bernoulli)
        
    z_mean = clamp(unsafeMean(q_z), tiny, 1.0 - tiny)

    outbound_dist.factors[1].m   =   unsafeMean(q_x)
    outbound_dist.factors[2].m   =   unsafeMean(q_x)
    outbound_dist.factors[1].W   =   z_mean*unsafeMean(q_w.factors[1])
    outbound_dist.factors[2].W   =   (1.0 - z_mean)*unsafeMean(q_w.factors[2])
    outbound_dist.factors[1].V   =   NaN
    outbound_dist.factors[2].V   =   NaN
    outbound_dist.factors[1].xi  =   NaN
    outbound_dist.factors[2].xi  =   NaN

    return outbound_dist
end
#
#
#
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

    z_mean = clamp(unsafeMean(q_z), tiny, 1.0 - tiny)

    outbound_dist.factors[1].m    =   unsafeMean(q_x)
    outbound_dist.factors[2].m    =   unsafeMean(q_x)
    outbound_dist.factors[1].W    =   z_mean*unsafeMean(q_w.factors[1])
    outbound_dist.factors[2].W    =   (1.0 - z_mean)*unsafeMean(q_w.factors[2])
    invalidate!(outbound_dist.factors[1].V)
    invalidate!(outbound_dist.factors[2].V)
    invalidate!(outbound_dist.factors[1].xi)
    invalidate!(outbound_dist.factors[2].xi)

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
        z_k_mean = clamp(z_mean[k], tiny, 1.0 - tiny)

        outbound_dist.factors[k].m = unsafeMean(q_x)
        outbound_dist.factors[k].W = z_k_mean*unsafeMean(q_w.factors[k])
        invalidate!(outbound_dist.factors[k].V)
        invalidate!(outbound_dist.factors[k].xi)
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
        z_k_mean = clamp(z_mean[k], tiny, 1.0 - tiny)

        outbound_dist.factors[k].m  = unsafeMean(q_x)
        outbound_dist.factors[k].W  = z_k_mean*unsafeMean(q_w.factors[k])
        outbound_dist.factors[k].V  = NaN
        outbound_dist.factors[k].xi = NaN
    end

    return outbound_dist
end
#
#
#
# VMP message towards i[:w]
# Univariate gaussian with two clusters
function variationalRule!{n_factors, T<:Univariate}(node::GaussianMixtureNode,
                                                    outbound_interface_index::Type{Val{2}},
                                                    outbound_dist::Partitioned{Gamma, n_factors},
                                                    q_m::Partitioned{T, n_factors},
                                                    q_w::Any,
                                                    q_x::Univariate,
                                                    q_z::Bernoulli)

    z_mean = clamp(unsafeMean(q_z), tiny, 1.0 - tiny)

    outbound_dist.factors[1].a   =   1.0 + 0.5*z_mean
    outbound_dist.factors[1].b   =   0.5*z_mean*( (unsafeMean(q_x) - unsafeMean(q_m.factors[1]))^2 + unsafeCov(q_x) + unsafeCov(q_m.factors[1]) )

    outbound_dist.factors[2].a   =   1.0 + 0.5*(1.0 - z_mean)
    outbound_dist.factors[2].b   =   0.5*(1.0 - z_mean)*( (unsafeMean(q_x) - unsafeMean(q_m.factors[2]))^2 + unsafeCov(q_x) + unsafeCov(q_m.factors[2]) )

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
                                                            outbound_dist::Partitioned{ForneyLab.Wishart{dims}, n_factors},
                                                            q_m::Partitioned{T, n_factors},
                                                            q_w::Any,
                                                            q_x::Multivariate{dims},
                                                            q_z::Bernoulli)

    z_mean = clamp(unsafeMean(q_z), tiny, 1.0 - tiny)

    outbound_dist.factors[1].nu = 1.0 + z_mean + dims
    outbound_dist.factors[1].V  = cholinv( z_mean*( (unsafeMean(q_x) - unsafeMean(q_m.factors[1]))*(unsafeMean(q_x) - unsafeMean(q_m.factors[1]))' + unsafeCov(q_x) + unsafeCov(q_m.factors[1]) ) )
    
    outbound_dist.factors[2].nu = 1.0 + (1.0 - z_mean) + dims
    outbound_dist.factors[2].V  = cholinv( (1.0 - z_mean)*( (unsafeMean(q_x) - unsafeMean(q_m.factors[2]))*(unsafeMean(q_x) - unsafeMean(q_m.factors[2]))' + unsafeCov(q_x) + unsafeCov(q_m.factors[2]) ) )

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

    # Calculate nu and V for each factor
    for k = 1:n_factors
        z_k_mean = clamp(z_mean[k], tiny, 1.0 - tiny)

        outbound_dist.factors[k].nu = 1.0 + z_k_mean + dims
        outbound_dist.factors[k].V  = cholinv( z_k_mean*( (unsafeMean(q_x) - unsafeMean(q_m.factors[k]))*(unsafeMean(q_x) - unsafeMean(q_m.factors[k]))' + unsafeCov(q_x) + unsafeCov(q_m.factors[k]) ) )
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

    # Calculate the values of a and b
    for k = 1:n_factors
        z_k_mean = clamp(z_mean[k], tiny, 1.0 - tiny)

        outbound_dist.factors[k].a = 1.0 + 0.5*z_k_mean
        outbound_dist.factors[k].b = 0.5*z_k_mean*( (unsafeMean(q_x) - unsafeMean(q_m.factors[k]))^2 + unsafeCov(q_x) + unsafeCov(q_m.factors[k]) )
    end

    return outbound_dist
end
#
#
#
# VMP message towards i[:z]
# Univariate gaussian with two clusters
function variationalRule!{n_factors, TN<:Univariate, TG<:Univariate}(   node::GaussianMixtureNode,
                                                                        outbound_interface_index::Type{Val{4}},
                                                                        outbound_dist::Bernoulli,
                                                                        q_m::Partitioned{TN, n_factors},
                                                                        q_w::Partitioned{TG, n_factors},
                                                                        q_x::Univariate,
                                                                        q_z::Any)

    rho = zeros(n_factors)
    for k = 1:n_factors
        rho[k] = clamp(exp(-averageEnergy(GaussianNode, q_m.factors[k], q_w.factors[k], q_x)), tiny, huge)
    end
    outbound_dist.p = rho[1]/sum(rho)

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

    rho = zeros(n_factors)
    for k = 1:n_factors
        rho[k] = clamp(exp(-averageEnergy(GaussianNode, q_m.factors[k], q_w.factors[k], q_x)), tiny, huge)
    end
    outbound_dist.p = rho[1]/sum(rho)

    return outbound_dist
end




# VMP message towards i[:z]
# Multivariate Gaussian with multiple clusters
function variationalRule!{dims,n_factors}(node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{4}},
                            outbound_dist::ForneyLab.Categorical{n_factors},
                            q_m::Partitioned{MvGaussian{dims},n_factors},
                            q_w::Partitioned{Wishart{dims}, n_factors},
                            q_x::MvGaussian{dims},
                            q_z::Any)

    ensureParameters!(q_x, (:m, :V))

    a = zeros(n_factors)
    ln_ro = zeros(n_factors)

    i = collect(1:dims)

    for k = 1:n_factors
      ensureParameters!(q_m.factors[k], (:m, :V))

      multidi = sum(digamma((q_w.factors[k].nu + 1 - i)/2))

      e_ln_w=  multidi + dims*log(2.0) + log(det(q_w.factors[k].V))

      e_w = q_w.factors[k].nu*q_w.factors[k].V

      gausterm = (transpose(q_x.m-q_m.factors[k].m)*e_w*(q_x.m - q_m.factors[k].m))[1] + trace((q_x.V + q_m.factors[k].V)*e_w)
      ln_ro[k]        =   0.5*e_ln_w - dims/2.0*log(2.0*pi) - 0.5*gausterm

    end

    sum_ro=sum(exp(ln_ro))

    if sum_ro > tiny
      for k = 1:n_factors
        outbound_dist.p[k] = exp(ln_ro[k])/sum_ro
      end
    else
      for k = 1:n_factors
        outbound_dist.p[k] = 1/n_factors
      end
    end

    return outbound_dist
end
#
# VMP message towards i[:z]
# Univariate gaussian with multiple clusters
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{4}},
                            outbound_dist::Categorical{n_factors},
                            q_m::Partitioned{Gaussian,n_factors},
                            q_w::Partitioned{Gamma,n_factors},
                            q_x::Gaussian,
                            q_z::Any)


      ensureParameters!(q_x, (:m, :V))

      a = zeros(n_factors)
      ln_ro = zeros(n_factors)


      #Calculate rho for each cluster
      for k = 1:n_factors
        ensureParameters!(q_m.factors[k], (:m, :V))

        e_ln_w     =   digamma(q_w.factors[k].a) - log(q_w.factors[k].b)
        e_m_square =   q_x.m^2 - 2.0*q_x.m*q_m.factors[k].m + q_m.factors[k].V + q_m.factors[k].m^2 + q_x.V
        ln_ro[k]   =   0.5*e_ln_w - 0.5*log(2pi) - 0.5*q_w.factors[k].a/q_w.factors[k].b*e_m_square

      end

      sum_ro=sum(exp(ln_ro))

      #normalize rho
      for k = 1:n_factors
          if sum_ro > tiny
            outbound_dist.p[k] = exp(ln_ro[k])/sum_ro
          else
            outbound_dist.p[k] = 1/n_factors
          end
      end



    return outbound_dist
end
#
#
# VMP message towards i[:x]
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::Gaussian,
                            q_m::Partitioned{Gaussian,n_factors},
                            q_w::Partitioned{Gamma,n_factors},
                            q_x::Any,
                            q_z::Bernoulli)



    V = pinv(q_z.p*q_w.factors[1].a/q_w.factors[1].b + (1. - q_z.p)*q_w.factors[2].a/q_w.factors[2].b)
    outbound_dist.m = (q_z.p*q_m.factors[1].m*q_w.factors[1].a/q_w.factors[1].b + (1. - q_z.p)*q_m.factors[2].m*q_w.factors[2].a/q_w.factors[2].b)*V
    outbound_dist.V = V
    outbound_dist.W = NaN
    outbound_dist.xi = NaN


    return outbound_dist
end
#
# #
# VMP message towards i[:x]
function variationalRule!{dims,n_factors}(  node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::MvGaussian{dims},
                            q_m::Partitioned{MvGaussian{dims}, n_factors},
                            q_w::Partitioned{ForneyLab.Wishart{dims},n_factors},
                            q_x::Any,
                            q_z::ForneyLab.Bernoulli)

    lambda = (q_z.p*q_w.factors[1].nu*q_w.factors[1].V + (1. - q_z.p)*q_w.factors[2].nu*q_w.factors[2].V)
    outbound_dist.xi = (q_z.p*q_w.factors[1].nu*q_w.factors[1].V*q_m.factors[1].m + (1. - q_z.p)*q_w.factors[2].nu*q_w.factors[2].V*q_m.factors[2].m)
    outbound_dist.W = (lambda)
    invalidate!(outbound_dist.m)
    invalidate!(outbound_dist.V)

    return outbound_dist
end
#
#
#
#
# #
# VMP message towards i[:x]
function variationalRule!{dims,n_factors}(  node::GaussianMixtureNode,
                            outbound_interface_index::Type{Val{3}},
                            outbound_dist::MvGaussian{dims},
                            q_m::Partitioned{MvGaussian{dims}, n_factors},
                            q_w::Partitioned{Wishart{dims},n_factors},
                            q_x::Any,
                            q_z::Categorical{n_factors})

    lambda = 0.0
    xi = 0.0

    for k = 1:n_factors
        lambda = lambda + q_z.p[k]*q_w.factors[k].nu*q_w.factors[k].V
        xi = xi + q_z.p[k]*q_w.factors[k].nu*q_w.factors[k].V*q_m.factors[k].m
    end

    outbound_dist.xi = xi
    outbound_dist.W = lambda
    invalidate!(outbound_dist.m)
    invalidate!(outbound_dist.V)

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

    lambda = 0.0
    m = 0.0

    for k = 1:n_factors
        lambda = lambda + q_z.p[k]*q_w.factors[k].a/q_w.factors[k].b
        m = m + q_z.p[k]*q_m.factors[k].m*q_w.factors[k].a/q_w.factors[k].b

    end

    outbound_dist.V = pinv(lambda)
    outbound_dist.m = m*outbound_dist.V
    outbound_dist.W = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end
