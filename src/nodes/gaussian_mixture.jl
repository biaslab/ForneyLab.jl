export GaussianMixtureNode

"""
Description:

    f(m1,m2,w1,w2,pi,x,z) = N(x|m1,1/w1)^z[1] * N(x|m2,1/w2)^z[2] * pi[1]^z[1] * pi[2]^z[2]

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
        self = new(id, Array(Interface, 7), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:pi, :m1, :w1, :m2, :w2, :x, :z])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::GaussianMixtureNode) = false


# VMP message towards i[:pi]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{1}},
                            outbound_dist::BetaDistribution,
                            ::Any,
                            q_m1::GaussianDistribution,
                            q_w1::GammaDistribution,
                            q_m2::GaussianDistribution,
                            q_w2::GammaDistribution,
                            q_x::GaussianDistribution,
                            q_z::BernoulliDistribution)

    ensureParameters!(q_m1, (:m, :V))


    outbound_dist.a = q_z.p
    outbound_dist.b = 1.-q_z.p

    return outbound_dist
end

# VMP message towards i[:m1]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{2}},
                            outbound_dist::GaussianDistribution,
                            q_pi::BetaDistribution,
                            ::Any,
                            q_w1::GammaDistribution,
                            q_m2::GaussianDistribution,
                            q_w2::GammaDistribution,
                            q_x::GaussianDistribution,
                            q_z::BernoulliDistribution)
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.m = q_x.m
    outbound_dist.V = NaN
    outbound_dist.xi  = NaN
    outbound_dist.W = q_z.p*q_w1.a/q_w1.b

    return outbound_dist
end

# VMP message towards i[:w1]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{3}},
                            outbound_dist::GammaDistribution,
                            q_pi::BetaDistribution,
                            q_m1::GaussianDistribution,
                            ::Any,
                            q_m2::GaussianDistribution,
                            q_w2::GammaDistribution,
                            q_x::GaussianDistribution,
                            q_z::BernoulliDistribution)

    ensureParameters!(q_m1, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.a = 1.+0.5*q_z.p
    e_m1_square = q_m1.V+q_m1.m^2
    outbound_dist.b = 0.5*q_z.p*(q_x.m^2-2.*q_x.m*q_m1.m+e_m1_square)

    return outbound_dist
end

# VMP message towards i[:m2]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{4}},
                            outbound_dist::GaussianDistribution,
                            q_pi::BetaDistribution,
                            q_m1::GaussianDistribution,
                            q_w1::GammaDistribution,
                            ::Any,
                            q_w2::GammaDistribution,
                            q_x::GaussianDistribution,
                            q_z::BernoulliDistribution)
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.m = q_x.m
    outbound_dist.V = NaN
    outbound_dist.xi  = NaN
    outbound_dist.W = (1.-q_z.p)*q_w2.a/q_w2.b

    return outbound_dist
end

# VMP message towards i[:w2]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{5}},
                            outbound_dist::GammaDistribution,
                            q_pi::BetaDistribution,
                            q_m1::GaussianDistribution,
                            q_w1::GammaDistribution,
                            q_m2::GaussianDistribution,
                            ::Any,
                            q_x::GaussianDistribution,
                            q_z::BernoulliDistribution)

    ensureParameters!(q_m2, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    outbound_dist.a = 1.+0.5*(1.-q_z.p)
    e_m2_square = q_m2.V+q_m2.m^2
    outbound_dist.b = 0.5*(1.-q_z.p)*(q_x.m^2-2.*q_x.m*q_m2.m+e_m2_square)

    return outbound_dist
end

# VMP message towards i[:z]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{7}},
                            outbound_dist::BernoulliDistribution,
                            q_pi::BetaDistribution,
                            q_m1::GaussianDistribution,
                            q_w1::GammaDistribution,
                            q_m2::GaussianDistribution,
                            q_w2::GammaDistribution,
                            q_x::GaussianDistribution,
                            ::Any)

    ensureParameters!(q_m1, (:m, :V))
    ensureParameters!(q_x, (:m, :V))

    e_ln_pi1 = digamma(q_pi.a)-digamma(q_pi.a+q_pi.b)
    e_ln_w1 = digamma(q_w1.a)-log(q_w1.b)
    e_m1_square = q_m1.V+q_m1.m^2
    ln_ro1 = e_ln_pi1+0.5*e_ln_w1-0.5*log(2pi)-0.5*q_w1.a/q_w1.b*e_m1_square
    outbound_dist.p = exp(ln_ro1)

    return outbound_dist
end

# VMP message towards i[:z]
function variationalRule!(  node::GaussianMixtureNode,
                            ::Type{Val{6}},
                            outbound_dist::GaussianDistribution,
                            q_pi::BetaDistribution,
                            q_m1::GaussianDistribution,
                            q_w1::GammaDistribution,
                            q_m2::GaussianDistribution,
                            q_w2::GammaDistribution,
                            q_x::Any,
                            q_z::BernoulliDistribution)

    outbound_dist.m=0.0
    outbound_dist.V=100
    outbound_dist.xi=NaN
    outbound_dist.W=NaN

    return outbound_dist
end
