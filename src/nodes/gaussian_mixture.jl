export GaussianMixtureNode

"""
Description:

    Univariate and multivariate gaussian mixture model for data with two clusters:
    f(m1,m2,w1,w2,pi,x,z) = N(x|m1,1/w1)^z[1] * N(x|m2,1/w2)^z[2] * pi[1]^z[1] * pi[2]^z[2]


                  | pi
                  |
          ________|________
      m1  |               |  x
     -----|               |----
      w1  |               |  z
     -----|       GM      |----
      m2  |               |
     -----|               |
      w2  |               |
     -----|               |
          |_______________|


Interfaces:

    1 i[:pi]
    2 i[:m1]
    3 i[:w1]
    4 i[:m2]
    5 i[:w2]
    6 i[:x]
    7 i[:z]

Construction:

    GaussianMixtureNode(id=:my_node)
"""
type GaussianMixtureNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function GaussianMixtureNode(; id=generateNodeId(GaussianMixtureNode))
        #self = new(id, Array(Interface, 5), Dict{Symbol,Interface}())
        self = new(id, Array(Interface, 7), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        #for (iface_index, iface_handle) in enumerate([:pi, :m, :w, :x, :z])
        for (iface_index, iface_handle) in enumerate([:pi, :m1, :w1, :m2, :w2, :x, :z])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::GaussianMixtureNode) = false



# VMP message towards i[:pi]
# Univariate gaussian with two clusters
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{1}},
                            outbound_dist::Beta,
                            ::Any,
                            q_m1::Gaussian,
                            q_w1::Gamma,
                            q_m2::Gaussian,
                            q_w2::Gamma,
                            q_x::Gaussian,
                            q_z::Bernoulli)

    ensureParameters!(q_m1, (:m, :V))


    outbound_dist.a   =   q_z.p+1.
    outbound_dist.b   =   2.-q_z.p

    return outbound_dist
end

# VMP message towards i[:pi]
# Multivariate gaussian with two clusters
function variationalRule!{dims}(  node::GaussianMixtureNode,
                            ::Type{Val{1}},
                            outbound_dist::Beta,
                            ::Any,
                            q_m1::MvGaussian{dims},
                            q_w1::Wishart{dims},
                            q_m2::MvGaussian{dims},
                            q_w2::Wishart{dims},
                            q_x::MvGaussian{dims},
                            q_z::Bernoulli)

    ensureParameters!(q_m1, (:m, :V))


     outbound_dist.a   =   q_z.p+1.
     outbound_dist.b   =   2.-q_z.p

    return outbound_dist
end


# VMP message towards i[:m1]
# Univariate Gaussian with two clusters
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{2}},
                            outbound_dist::Gaussian,
                            q_pi::Beta,
                            ::Any,
                            q_w1::Gamma,
                            q_m2::Gaussian,
                            q_w2::Gamma,
                            q_x::Gaussian,
                            q_z::Bernoulli)
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.m   =   q_x.m
    outbound_dist.V   =   NaN
    outbound_dist.xi  =   NaN
    outbound_dist.W   =   q_z.p*q_w1.a/q_w1.b

    return outbound_dist
end

# VMP message towards i[:m1]
# Multivariate Gaussian with two clusters
function variationalRule!{dims}(  node::GaussianMixtureNode,
                            ::Type{Val{2}},
                            outbound_dist::MvGaussian{dims},
                            q_pi::Beta,
                            ::Any,
                            q_w1::Wishart{dims},
                            q_m2::MvGaussian{dims},
                            q_w2::Wishart{dims},
                            q_x::MvGaussian{dims},
                            q_z::Bernoulli)
    ensureParameters!(q_x, (:m, :V))

     outbound_dist.m    =   deepcopy(q_x.m)
     invalidate!(outbound_dist.V)
     invalidate!(outbound_dist.xi)
     outbound_dist.W    =   q_z.p*q_w1.nu*q_w1.V

    return outbound_dist
end


# VMP message towards i[:w1]
# Univariate Gaussian with two clusters
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{3}},
                            outbound_dist::Gamma,
                            q_pi::Beta,
                            q_m1::Gaussian,
                            ::Any,
                            q_m2::Gaussian,
                            q_w2::Gamma,
                            q_x::Gaussian,
                            q_z::Bernoulli)

    ensureParameters!(q_m1, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.a   =   1.+0.5*q_z.p
    e_m1_square       =   q_m1.V+q_m1.m^2
    outbound_dist.b   =   0.5*q_z.p*(q_x.m^2-2.*q_x.m*q_m1.m+e_m1_square)

    return outbound_dist
end

# VMP message towards i[:w1]
# Multivariate Gaussian with two clusters
function variationalRule!{dims}(  node::GaussianMixtureNode,
                            ::Type{Val{3}},
                            outbound_dist::Wishart{dims},
                            q_pi::Beta,
                            q_m1::MvGaussian{dims},
                            ::Any,
                            q_m2::MvGaussian{dims},
                            q_w2::Wishart{dims},
                            q_x::MvGaussian{dims},
                            q_z::Bernoulli)

    ensureParameters!(q_m1, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.nu    =   1.+q_z.p+dims
    gausterm1           =   (deepcopy(q_x.m)-q_m1.m)*transpose(deepcopy(q_x.m)-q_m1.m)+q_m1.V

    #if statement to prevent multiplication with zero
    if det((q_z.p)*gausterm1)<tiny
      outbound_dist.V     =   pinv((q_z.p)*gausterm1+eye(dims)*tiny)
    else
      outbound_dist.V     =   pinv(q_z.p*gausterm1)
    end

    return outbound_dist
end


# VMP message towards i[:m2]
# Univariate Gaussian with two clusters
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{4}},
                            outbound_dist::Gaussian,
                            q_pi::Beta,
                            q_m1::Gaussian,
                            q_w1::Gamma,
                            ::Any,
                            q_w2::Gamma,
                            q_x::Gaussian,
                            q_z::Bernoulli)
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.m   =   q_x.m
    outbound_dist.V   =   NaN
    outbound_dist.xi  =   NaN
    outbound_dist.W   =   (1.-q_z.p)*q_w2.a/q_w2.b

    return outbound_dist
end

# VMP message towards i[:m2]
# Multivariate Gaussian with two clusters
function variationalRule!{dims}(  node::GaussianMixtureNode,
                            ::Type{Val{4}},
                            outbound_dist::MvGaussian{dims},
                            q_pi::Beta,
                            q_m1::MvGaussian{dims},
                            q_w1::Wishart{dims},
                            ::Any,
                            q_w2::Wishart{dims},
                            q_x::MvGaussian{dims},
                            q_z::Bernoulli)
     ensureParameters!(q_x, (:m, :V))

     outbound_dist.m   =   deepcopy(q_x.m)

     invalidate!(outbound_dist.V)
     invalidate!(outbound_dist.xi)
    outbound_dist.W   =   (1.-q_z.p)*q_w2.nu*q_w2.V

    return outbound_dist
end

# VMP message towards i[:w2]
# Univariate Gaussian with two clusters
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{5}},
                            outbound_dist::Gamma,
                            q_pi::Beta,
                            q_m1::Gaussian,
                            q_w1::Gamma,
                            q_m2::Gaussian,
                            ::Any,
                            q_x::Gaussian,
                            q_z::Bernoulli)

    ensureParameters!(q_m2, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.a   =   1.+0.5*(1.-q_z.p)
    e_m2_square       =   q_m2.V+q_m2.m^2
    outbound_dist.b   =   0.5*(1.-q_z.p)*(q_x.m^2-2.*q_x.m*q_m2.m+e_m2_square)
    return outbound_dist
end
#
# VMP message towards i[:w2]
# Multivariate Gaussian with two clusters
function variationalRule!{dims}(  node::GaussianMixtureNode,
                            ::Type{Val{5}},
                            outbound_dist::Wishart{dims},
                            q_pi::Beta,
                            q_m1::MvGaussian{dims},
                            q_w1::Wishart{dims},
                            q_m2::MvGaussian{dims},
                            ::Any,
                            q_x::MvGaussian{dims},
                            q_z::Bernoulli)

    ensureParameters!(q_m2, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.nu    =   1.+(1.-q_z.p)+dims
    gausterm2           =   (deepcopy(q_x.m)-q_m2.m)*transpose(deepcopy(q_x.m)-q_m2.m)+q_m2.V

    #if statement to prevent multiplication by zero

    if det((1.-q_z.p)*gausterm2)<tiny
      outbound_dist.V     =   pinv((1.-q_z.p)*gausterm2+eye(dims)*tiny)

    else
      outbound_dist.V     =   pinv((1.-q_z.p)*gausterm2)

    end

    return outbound_dist
end


# VMP message towards i[:z]
# Univariate Gaussian with two clusters
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{7}},
                            outbound_dist::Bernoulli,
                            q_pi::Beta,
                            q_m1::Gaussian,
                            q_w1::Gamma,
                            q_m2::Gaussian,
                            q_w2::Gamma,
                            q_x::Gaussian,
                            ::Any)

    ensureParameters!(q_m1, (:m, :V))
    ensureParameters!(q_x, (:m, :V))
    ensureParameters!(q_m2, (:m, :V))

    #calculating ln(ro1)
    e_ln_pi1    =   digamma(q_pi.a)-digamma(q_pi.a+q_pi.b)
    e_ln_w1     =   digamma(q_w1.a)-log(q_w1.b)
    e_m1_square =   q_x.m^2-2.0*q_x.m*q_m1.m+q_m1.V+q_m1.m^2
    ln_ro1      =   e_ln_pi1+0.5*e_ln_w1-0.5*log(2pi)-0.5*q_w1.a/q_w1.b*e_m1_square

    #calculating ln(ro2) for normalization

    e_ln_pi2    =   digamma(q_pi.b)-digamma(q_pi.b+q_pi.a)
    e_ln_w2     =   digamma(q_w2.a)-log(q_w2.b)
    e_m2_square =   q_x.m^2-2.0*q_x.m*q_m2.m+q_m2.V+q_m2.m^2
    ln_ro2      =   e_ln_pi2+0.5*e_ln_w2-0.5*log(2pi)-0.5*q_w2.a/q_w2.b*e_m2_square

    outbound_dist.p = exp(ln_ro1)/(exp(ln_ro1)+exp(ln_ro2))

    return outbound_dist
end

# VMP message towards i[:z]
# Multivariate Gaussian with two clusters
function variationalRule!{dims}(node::GaussianMixtureNode,
                            ::Type{Val{7}},
                            outbound_dist::Bernoulli,
                            q_pi::Beta,
                            q_m1::MvGaussian{dims},
                            q_w1::Wishart{dims},
                            q_m2::MvGaussian{dims},
                            q_w2::Wishart{dims},
                            q_x::MvGaussian{dims},
                            ::Any)

    ensureParameters!(q_m1, (:m, :V))
    ensureParameters!(q_x, (:m, :V))
    ensureParameters!(q_m2, (:m, :V))

    e_ln_pi1      =   digamma(q_pi.a)-digamma(q_pi.a+q_pi.b)

    #multivariate digamma
    i = collect(1:dims)
    multidi1=sum(digamma((q_w1.nu+1-i)/2))

    e_ln_w1       =   deepcopy(multidi1)+dims*log(2.0)+log(det(q_w1.V))
    e_w1          =   q_w1.nu*q_w1.V
    gausterm1 = (transpose(q_x.m-q_m1.m)*e_w1*(q_x.m-q_m1.m))[1] + trace(q_m1.V*e_w1)

    ln_ro1        =   e_ln_pi1+0.5*e_ln_w1-dims/2.0*log(2.0*pi)-0.5*gausterm1


    #calculating ln(ro2) for normalization

    e_ln_pi2      =   digamma(q_pi.b)-digamma(q_pi.b+q_pi.a)

    #multivariate digamma
    i = collect(1:dims)
    multidi2= sum(digamma((q_w2.nu+1-i)/2))

    e_ln_w2       =  multidi2+dims*log(2.0)+log(det(q_w2.V))
    e_w2          =   q_w2.nu*q_w2.V
    gausterm2     =  (transpose(q_x.m-q_m2.m)*e_w2*(q_x.m-q_m2.m))[1] + trace(q_m2.V*e_w2) # wordt te groot


    ln_ro2        =   e_ln_pi2+0.5*e_ln_w2-dims/2.0*log(2.0*pi)-0.5*gausterm2

    #Normalize message
    #if statement to prevent division by zero
    if exp(ln_ro1)+exp(ln_ro2)>tiny
      outbound_dist.p = exp(ln_ro1)/(exp(ln_ro1)+exp(ln_ro2))
    else
      outbound_dist.p=0.5
    end
    return outbound_dist
end




# VMP message towards i[:x]
function variationalRule!(  node::GaussianMixtureNode,s
                            ::Type{Val{6}},
                            outbound_dist::Gaussian,
                            q_pi::Beta,
                            q_m1::Gaussian,
                            q_w1::Gamma,
                            q_m2::Gaussian,
                            q_w2::Gamma,
                            q_x::Any,
                            q_z::Bernoulli)

    outbound_dist.m=0.0
    outbound_dist.V=100
    outbound_dist.xi=NaN
    outbound_dist.W=NaN

    return outbound_dist
end

# VMP message towards i[:x]
function variationalRule!{dims}(  node::GaussianMixtureNode,
                            ::Type{Val{6}},
                            outbound_dist::MvGaussian{dims},
                            q_pi::Beta,
                            q_m1::MvGaussian{dims},
                            q_w1::Wishart{dims},
                            q_m2::MvGaussian{dims},
                            q_w2::Wishart{dims},
                            q_x::Any,
                            q_z::Bernoulli)

    outbound_dist.m   =  ones(dims)
    outbound_dist.V   = 100.0*eye(dims)
    invalidate!(outbound_dist.xi)
    invalidate!(outbound_dist.W)

    return outbound_dist
end
