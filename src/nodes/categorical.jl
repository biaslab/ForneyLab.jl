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
                            pi::Dirichlet{n_factors},
                            z::Any)

    k = collect(1:n_factors)
    sum_a = sum(q_pi.alpha[k])
    ro = zeros(n_factors)

    for k = 1:n_factors
      e_ln_pi   =   digamma(q_pi.alpha[k]) - digamma(sum_a)
      ro[k]   =   exp(e_ln_pi)
    end

    sum_ro=sum(ro)

    #normalize rho
    for k = 1:n_factors
      if sum_ro > tiny
        outbound_dist.p[k] = ro[k]/sum_ro
      else
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
function variationalRule!(  node::CategoricalNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Dirichlet{n_factors},
                            pi::Any,
                            z::Categorical{n_factors})


    outbound_dist.alpha = z.p+1.

    return outbound_dist
end
