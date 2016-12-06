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
#                             msg_m::Message{PartitionedDistribution{Gaussian,2}},
#                             msg_w::Message{PartitionedDistribution{Gamma,2}},
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
#                             msg_m::Message{PartitionedDistribution{Gaussian,n_factors}},
#                             msg_w::Message{PartitionedDistribution{Gamma,n_factors}},
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
#                             msg_m::Message{PartitionedDistribution{MvGaussian{dims},n_factors}},
#                             msg_w::Message{PartitionedDistribution{Wishart{dims},n_factors}},
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
#                             msg_m::Message{PartitionedDistribution{MvGaussian{dims},n_factors}},
#                             msg_w::Message{PartitionedDistribution{Wishart{dims},n_factors}},
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
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{1}},
                            outbound_dist::Partitioned{Gaussian,n_factors},
                            ::Any,
                            q_w::Partitioned{Gamma,n_factors},
                            q_x::Gaussian,
                            q_z::Bernoulli)
        ensureParameters!(q_x, (:m,:V))

        outbound_dist.factors[1].m   =   q_x.m
        outbound_dist.factors[2].m   =   q_x.m
        outbound_dist.factors[1].V   =   NaN
        outbound_dist.factors[1].xi  =   NaN
        outbound_dist.factors[2].V   =   NaN
        outbound_dist.factors[2].xi  =   NaN
        outbound_dist.factors[1].W   =   q_z.p*q_w.factors[1].a/q_w.factors[1].b
        outbound_dist.factors[2].W   =   (1. - q_z.p)*q_w.factors[2].a/q_w.factors[2].b

    return outbound_dist
end
#
#
#
# VMP message towards i[:m]
# Multivariate Gaussian with two clusters
function variationalRule!{dims,n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{1}},
                            outbound_dist::Partitioned{MvGaussian{dims},n_factors},
                            ::Any,
                            q_w::Partitioned{ForneyLab.Wishart{dims},n_factors},
                            q_x::MvGaussian{dims},
                            q_z::ForneyLab.Bernoulli)

    ensureParameters!(q_x, (:m, :V))

     outbound_dist.factors[1].m    =   deepcopy(q_x.m)
     outbound_dist.factors[2].m    =   deepcopy(q_x.m)
     invalidate!(outbound_dist.factors[1].V)
     invalidate!(outbound_dist.factors[2].V)
     invalidate!(outbound_dist.factors[1].xi)
     invalidate!(outbound_dist.factors[2].xi)
     outbound_dist.factors[1].W    =   (q_z.p*q_w.factors[1].nu*q_w.factors[1].V)
     outbound_dist.factors[2].W    =   ((1. - q_z.p)*q_w.factors[2].nu*q_w.factors[2].V)

    return outbound_dist
end
#
#
#
#
# VMP message towards i[:m]
# Multivariate Gaussian with multiple clusters
function variationalRule!{dims,n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{1}},
                            outbound_dist::Partitioned{MvGaussian{dims},n_factors},
                            ::Any,
                            q_w::Partitioned{Wishart{dims},n_factors},
                            q_x::MvGaussian{dims},
                            q_z::Categorical{n_factors})

    ensureParameters!(q_x, (:m, :V))


    for k = 1:n_factors
      outbound_dist.factors[k].m  = deepcopy(q_x.m)
      invalidate!(outbound_dist.factors[k].V)
      invalidate!(outbound_dist.factors[k].xi)
      outbound_dist.factors[k].W  = q_z.p[k]*q_w.factors[k].nu*q_w.factors[k].V
    end

    return outbound_dist
end
#
# # VMP message towards i[:m]
# # Univariate Gaussian with multiple clusters
# function variationalRule!{n_factors}(  node::GaussianMixtureNode,
#                             ::Type{Val{2}},
#                             outbound_dist::Partitioned{Gaussian,n_factors},
#                             q_pi::Dirichlet{n_factors},
#                             ::Any,
#                             q_w::Partitioned{Gamma,n_factors},
#                             q_x::Gaussian,
#                             q_z::Categorical{n_factors})
#
#     ensureParameters!(q_x, (:m, :V))
#
#     for k=1:n_factors
#        outbound_dist.factors[k].m  = deepcopy(q_x.m)
#        outbound_dist.factors[k].V  = NaN
#        outbound_dist.factors[k].xi = NaN
#        outbound_dist.factors[k].W  = q_z.p[k]*q_w.factors[k].a/q_w.factors[k].b
#     end
#
#     return outbound_dist
# end
#
#
#
# VMP message towards i[:w]
# Univariate gaussian with two clusters
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{2}},
                            outbound_dist::Partitioned{Gamma,n_factors},
                            q_m::Partitioned{Gaussian,n_factors},
                            ::Any,
                            q_x::Gaussian,
                            q_z::Bernoulli)

      ensureParameters!(q_m.factors[1], (:m, :V))
      ensureParameters!(q_m.factors[2], (:m, :V))
      ensureParameters!(q_x, (:m, :V))

      outbound_dist.factors[1].a   =   1. + 0.5*q_z.p
      e_m1_square                  =   q_m.factors[1].V + q_m.factors[1].m^2 + q_x.V
      outbound_dist.factors[1].b   =   0.5*q_z.p*(q_x.m^2-2.*q_x.m*q_m.factors[1].m + e_m1_square)

      outbound_dist.factors[2].a   =   1. + 0.5*(1. - q_z.p)
      e_m2_square                  =   q_m.factors[2].V + q_m.factors[2].m^2 + q_x.V
      outbound_dist.factors[2].b   =   0.5*(1. - q_z.p)*(q_x.m^2 - 2.*q_x.m*q_m.factors[2].m + e_m2_square)

    return outbound_dist
end
#
#
#
# VMP message towards i[:w]
# Multivariate Gaussian with two clusters
function variationalRule!{dims,n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{2}},
                            outbound_dist::Partitioned{ForneyLab.Wishart{dims},n_factors},
                            q_m::Partitioned{MvGaussian{dims},n_factors},
                            ::Any,
                            q_x::MvGaussian{dims},
                            q_z::ForneyLab.Bernoulli)

    #Ensure that the messages have the right parameters
    ensureParameters!(q_m.factors[1], (:m, :V))
    ensureParameters!(q_x, (:m, :V))
    ensureParameters!(q_m.factors[2], (:m, :V))

    #Calculate nu and C for the first factor
    outbound_dist.factors[1].nu = 1. + q_z.p + dims
    gausterm1  =  (deepcopy(q_x.m) - q_m.factors[1].m)*transpose(deepcopy(q_x.m) - q_m.factors[1].m) + (q_m.factors[1].V) + (q_x.V)

    #if statement to prevent multiplication with zero
    if det((q_z.p)*gausterm1) < tiny
      outbound_dist.factors[1].V = pinv((q_z.p)*gausterm1 + eye(dims)*tiny)
    else
      outbound_dist.factors[1].V = pinv(q_z.p*gausterm1)
    end



    outbound_dist.factors[2].nu  = 1. + (1. - q_z.p) + dims
    gausterm2 = (deepcopy(q_x.m) - q_m.factors[2].m)*transpose(deepcopy(q_x.m) - q_m.factors[2].m) + q_m.factors[2].V + q_x.V

    #if statement to prevent multiplication by zero

    if det((1. - q_z.p)*gausterm2) < tiny
      outbound_dist.factors[2].V     =   pinv((1. - q_z.p)*gausterm2 + eye(dims)*tiny)

    else
      outbound_dist.factors[2].V     =   pinv((1. - q_z.p)*gausterm2)

    end
    return outbound_dist
end

#
#
# VMP message towards i[:w]
# Multivariate Gaussian with multiple clusters
function variationalRule!{dims,n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{2}},
                            outbound_dist::Partitioned{Wishart{dims},n_factors},
                            q_m::Partitioned{MvGaussian{dims},n_factors},
                            ::Any,
                            q_x::MvGaussian{dims},
                            q_z::Categorical{n_factors})

    #Ensure that the distributions have the correct parameters
    ensureParameters!(q_x,(:m,:V))

    #Calculate nu and V for each distribution
    for k = 1:n_factors
      ensureParameters!(q_m.factors[k], (:m, :V))

      outbound_dist.factors[k].nu = 1. + q_z.p[k] + dims
      gausterm = (deepcopy(q_x.m) - q_m.factors[k].m)*transpose(deepcopy(q_x.m) - q_m.factors[k].m) + q_m.factors[k].V + q_x.V

      #if statement to prevent multiplication with zero
      if det((q_z.p[k])*gausterm) < tiny
        outbound_dist.factors[k].V = pinv((q_z.p[k])*gausterm + eye(dims)*tiny)
      else
        outbound_dist.factors[k].V = pinv(q_z.p[k]*gausterm)
      end
    end

    return outbound_dist
end

# # VMP message towards i[:w]
# # Univariate Gaussian with multiple clusters
# function variationalRule!{n_factors}(  node::GaussianMixtureNode,
#                             ::Type{Val{3}},
#                             outbound_dist::Partitioned{Gamma,n_factors},
#                             q_pi::Dirichlet{n_factors},
#                             q_m::Partitioned{Gaussian,n_factors},
#                             ::Any,
#                             q_x::Gaussian,
#                             q_z::Categorical{n_factors})
#
#     #Ensure that the distribution has the correct parameters
#     ensureParameters!(q_x, (:m, :V))
#
#     # Calculate the values of a and b
#     for k = 1:n_factors
#       ensureParameters!(q_m.factors[k], (:m, :V))
#
#       outbound_dist.factors[k].a  = 1. + 0.5*q_z.p[k]
#       e_m1_square                 = q_m.factors[k].V + q_m.factors[k].m^2 + q_x.V
#       outbound_dist.factors[k].b  = 0.5*q_z.p[k]*(q_x.m^2 - 2.*q_x.m*q_m.factors[k].m + e_m1_square)
#
#       if outbound_dist.factors[k].b < tiny
#         outbound_dist.factors[k].b = 100*tiny + 0.5*q_z.p[k]*(q_x.m^2 - 2.*q_x.m*q_m.factors[k].m + e_m1_square)
#       end
#
#     end
#
#
#     return outbound_dist
# end
#
#
#
# VMP message towards i[:z]
# Univariate gaussian with two clusters
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{4}},
                            outbound_dist::Bernoulli,
                            q_m::Partitioned{Gaussian,n_factors},
                            q_w::Partitioned{Gamma,n_factors},
                            q_x::Gaussian,
                            ::Any)

      ensureParameters!(q_m.factors[1], (:m, :V))
      ensureParameters!(q_x, (:m, :V))
      ensureParameters!(q_m.factors[2], (:m, :V))

      #calculating ln(ro1)
      #e_ln_pi1    =   digamma(q_pi.a) - digamma(q_pi.a + q_pi.b)
      e_ln_w1     =   digamma(q_w.factors[1].a) - log(q_w.factors[1].b)
      e_m1_square =   q_x.m^2 - 2.0*q_x.m*q_m.factors[1].m + q_m.factors[1].V + q_m.factors[1].m^2 + q_x.V
      ln_ro1      =   0.5*e_ln_w1 - 0.5*log(2pi) - 0.5*q_w.factors[1].a/q_w.factors[1].b*e_m1_square

      #calculating ln(ro2) for normalization

      #e_ln_pi2    =   digamma(q_pi.b) - digamma(q_pi.b + q_pi.a)
      e_ln_w2     =   digamma(q_w.factors[2].a) - log(q_w.factors[2].b)
      e_m2_square =   q_x.m^2 - 2.0*q_x.m*q_m.factors[2].m + q_m.factors[2].V + q_m.factors[2].m^2 + q_x.V
      ln_ro2      =   0.5*e_ln_w2 - 0.5*log(2pi) - 0.5*q_w.factors[2].a/q_w.factors[2].b*e_m2_square

      if (exp(ln_ro1) + exp(ln_ro2)) > tiny
          outbound_dist.p = exp(ln_ro1)/(exp(ln_ro1) + exp(ln_ro2))
      else
          outbound_dist.p = 1/2
      end


    return outbound_dist
end


# VMP message towards i[:z]
# Multivariate Gaussian with two clusters
function variationalRule!{dims,n_factors}(node::GaussianMixtureNode,
                            ::Type{Val{4}},
                            outbound_dist::ForneyLab.Bernoulli,
                            q_m::Partitioned{MvGaussian{dims},n_factors},
                            q_w::Partitioned{ForneyLab.Wishart{dims}, n_factors},
                            q_x::MvGaussian{dims},
                            ::Any)

    ensureParameters!(q_m.factors[1], (:m, :V))
    ensureParameters!(q_x, (:m, :V))
    ensureParameters!(q_m.factors[2], (:m, :V))

    # e_ln_pi1      =   digamma(q_pi.a) - digamma(q_pi.a + q_pi.b)

    #multivariate digamma
    i = collect(1:dims)
    multidi1=sum(digamma((q_w.factors[1].nu + 1 - i)/2))

    e_ln_w1       =   deepcopy(multidi1) + dims*log(2.0) + log(det(q_w.factors[1].V))
    e_w1          =   q_w.factors[1].nu*q_w.factors[1].V
    gausterm1 = (transpose(q_x.m - q_m.factors[1].m)*e_w1*(q_x.m - q_m.factors[1].m))[1] + trace((q_x.V + q_m.factors[1].V)*e_w1)

    ln_ro1        =   0.5*e_ln_w1-dims/2.0*log(2.0*pi) - 0.5*gausterm1


    #calculating ln(ro2) for normalization
    #e_ln_pi2      =   digamma(q_pi.b) - digamma(q_pi.b + q_pi.a)

    #multivariate digamma
    i = collect(1:dims)
    multidi2= sum(digamma((q_w.factors[2].nu + 1 - i)/2))

    e_ln_w2       =  multidi2 + dims*log(2.0) + log(det(q_w.factors[2].V))
    e_w2          =   q_w.factors[2].nu*q_w.factors[2].V
    gausterm2     =  (transpose(q_x.m - q_m.factors[2].m)*e_w2*(q_x.m - q_m.factors[2].m))[1] + trace((q_x.V + q_m.factors[2].V)*e_w2)


    ln_ro2        =   0.5*e_ln_w2 - dims/2.0*log(2.0*pi) - 0.5*gausterm2

    #Normalize message
    #if statement to prevent division by zero
    if exp(ln_ro1) + exp(ln_ro2) > tiny
      outbound_dist.p = exp(ln_ro1)/(exp(ln_ro1) + exp(ln_ro2))
    else
      outbound_dist.p = 0.5
    end
    return outbound_dist
end




# VMP message towards i[:z]
# Multivariate Gaussian with multiple clusters
function variationalRule!{dims,n_factors}(node::GaussianMixtureNode,
                            ::Type{Val{4}},
                            outbound_dist::Categorical{n_factors},
                            q_m::Partitioned{MvGaussian{dims},n_factors},
                            q_w::Partitioned{Wishart{dims}, n_factors},
                            q_x::MvGaussian{dims},
                            ::Any)

    ensureParameters!(q_x, (:m, :V))

    a = zeros(n_factors)
    ln_ro = zeros(n_factors)

    k = collect(1:n_factors)
    sum_a = sum(q_pi.alpha[k])

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
# # VMP message towards i[:z]
# # Univariate gaussian with multiple clusters
# function variationalRule!{n_factors}(  node::GaussianMixtureNode,
#                             ::Type{Val{5}},
#                             outbound_dist::Categorical{n_factors},
#                             q_pi::Dirichlet{n_factors},
#                             q_m::Partitioned{Gaussian,n_factors},
#                             q_w::Partitioned{Gamma,n_factors},
#                             q_x::Gaussian,
#                             ::Any)
#
#
#       ensureParameters!(q_x, (:m, :V))
#
#       a = zeros(n_factors)
#       ln_ro = zeros(n_factors)
#
#       k = collect(1:n_factors)
#       sum_a = sum(q_pi.alpha[k])
#
#
#       #Calculate rho for each cluster
#       for k = 1:n_factors
#         ensureParameters!(q_m.factors[k], (:m, :V))
#
#         e_ln_pi    =   digamma(q_pi.alpha[k]) - digamma(sum_a)
#         e_ln_w     =   digamma(q_w.factors[k].a) - log(q_w.factors[k].b)
#         e_m_square =   q_x.m^2 - 2.0*q_x.m*q_m.factors[k].m + q_m.factors[k].V + q_m.factors[k].m^2 + q_x.V
#         ln_ro[k]   =   e_ln_pi + 0.5*e_ln_w - 0.5*log(2pi) - 0.5*q_w.factors[k].a/q_w.factors[k].b*e_m_square
#
#       end
#
#       sum_ro=sum(exp(ln_ro))
#
#       #normalize rho
#       for k = 1:n_factors
#           if sum_ro > tiny
#             outbound_dist.p[k] = exp(ln_ro[k])/sum_ro
#           else
#             outbound_dist.p[k] = 1/n_factors
#           end
#       end
#
#
#
#     return outbound_dist
# end
#
#
# VMP message towards i[:x]
function variationalRule!{n_factors}(  node::GaussianMixtureNode,
                            ::Type{Val{3}},
                            outbound_dist::Gaussian,
                            q_m::Partitioned{Gaussian,n_factors},
                            q_w::Partitioned{Gamma,n_factors},
                            ::Any,
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
                            ::Type{Val{3}},
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
                            ::Type{Val{3}},
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
# # VMP message towards i[:x]
# function variationalRule!{n_factors}(  node::GaussianMixtureNode,
#                             ::Type{Val{4}},
#                             outbound_dist::Gaussian,
#                             q_pi::Dirichlet{n_factors},
#                             q_m::PartitionedDistribution{Gaussian,n_factors},
#                             q_w::PartitionedDistribution{Gamma,n_factors},
#                             ::Any,
#                             q_z::Categorical{n_factors})
#
#     lambda = 0.0
#     m = 0.0
#
#     for k = 1:n_factors
#         lambda = lambda + q_z.p[k]*q_w.factors[k].a/q_w.factors[k].b
#         m = m + q_z.p[k]*q_m.factors[k].m*q_w.factors[k].a/q_w.factors[k].b
#
#     end
#
#     outbound_dist.V = pinv(lambda)
#     outbound_dist.m = m*outbound_dist.V
#     outbound_dist.W = NaN
#     outbound_dist.xi = NaN
#
#     return outbound_dist
# end
