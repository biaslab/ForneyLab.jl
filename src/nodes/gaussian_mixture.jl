export GaussianMixtureNode

"""
Description:

    f(m1,m2,w1,w2,pi,x,z) = N(x|m1,1/w1)^z[1] * N(x|m2,1/w2)^z[2] * pi[1] * p[2]

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
                            q_x::DeltaDistribution{Float64},
                            q_z::BernoulliDistribution)

    ensureParameters!(q_m1, (:m, :V))

    q_m1.m + q_m1.V

    outbound_dist.a = 0.0
    outbound_dist.b = 0.0

    return outbound_dist
end
