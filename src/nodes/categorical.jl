export CategoricalNode

"""
Description:

    Categorical factor node: Cat(z|pi)
    Derivations can be found in the derivations document.

    in, pi       out, z
    ----->[Cat]]----->

    z ∈ {0, 1}
    pi ∈ [0, 1]

    f(z,pi) = pi_1^z_1 pi_2^z_2

Interfaces:

    1 i[:pi], 2 i[:z]

Construction:

    CategoricalNode(id=:my_node)
"""
type CategoricalNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function CategoricalNode(; id=generateNodeId(CategoricalNode))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:pi, :z])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::CategoricalNode) = false


############################################
# Variational update functions
############################################

"""
CategoricalNode:

    Dir       Cat
    --->[Cat]--->
              -->
"""
function variationalRule!{n_factors}(  node::CategoricalNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Categorical{n_factors},
                            q_pi::Dirichlet{n_factors},
                            q_z::Any)


    ro  =   exp(digamma(q_pi.alpha))

    sum_ro = sum(ro)

    #normalize rho
    if sum_ro > tiny
      outbound_dist.p= ro/sum_ro
    else
        for k = 1:n_factors
          outbound_dist.p[k] = 1/n_factors
        end
    end

    return outbound_dist
end

"""
CategoricalNode:

    Dir      Cat
    --->[Cat]--->
    <--
"""
function variationalRule!{n_factors}(  node::CategoricalNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Dirichlet{n_factors},
                            q_pi::Any,
                            q_z::Categorical{n_factors})


    outbound_dist.alpha = q_z.p+1.

    return outbound_dist
end

############################
# Average energy functional
############################

function averageEnergy(::Type{CategoricalNode}, q_pi::Dirichlet, q_z::Categorical)
    digamma(sum(q_pi.alpha))-(q_z.p'*digamma(q_pi.alpha))[1]
end
