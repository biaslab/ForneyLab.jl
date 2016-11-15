export BernoulliNode

"""
Description:

    Bernoulli factor node: Bern(x|c)
    Derivations can be found in the derivations document.

    in, x       out, c
    ----->[Bern]----->

    c ∈ {0, 1}
    x ∈ [0, 1]

    f(c,x) = x^c (1 - x)^(1 - c)

Interfaces:

    1 i[:in], 2 i[:out]

Construction:

    BernoulliNode(id=:my_node)
"""
type BernoulliNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function BernoulliNode(; id=generateNodeId(BernoulliNode))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:in, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::BernoulliNode) = false


############################################
# Sum-product update functions
############################################

"""
BernoulliNode:

    Beta      Bern
    --->[Bern]--->
              -->
"""
function sumProductRule!(   node::BernoulliNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Bernoulli,
                            msg_in::Message{Beta},
                            msg_out::Any)

    outbound_dist.p = msg_in.payload.a/(msg_in.payload.a + msg_in.payload.b)

    return outbound_dist
end

"""
BernoulliNode:

      δ       Bern
    --->[Bern]--->
              -->
"""
function sumProductRule!(   node::BernoulliNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Bernoulli,
                            msg_in::Message{Delta{Float64}},
                            msg_out::Any)

    outbound_dist.p = msg_in.payload.m

    return outbound_dist
end

"""
BernoulliNode:

    Beta       δ
    --->[Bern]--->
    <--
"""
function sumProductRule!(   node::BernoulliNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Beta,
                            msg_in::Any,
                            msg_out::Message{Delta{Bool}})

    outbound_dist.a = msg_out.payload.m + 1.0
    outbound_dist.b = 2.0 - msg_out.payload.m

    return outbound_dist
end


############################################
# Variational update functions
############################################

"""
BernoulliNode:

    Beta      Bern
    --->[Bern]--->
              -->
"""
function variationalRule!(  node::BernoulliNode,
                            outbound_interface_index::Type{Val{2}},
                            outbound_dist::Bernoulli,
                            marg_in::Beta,
                            marg_out::Any)

    if marg_in.a == marg_in.b # Catch for numeric reasons
        outbound_dist.p = 0.5
    else
        rho_a = exp(digamma(marg_in.a))
        rho_b = exp(digamma(marg_in.b))

        outbound_dist.p = rho_a/(rho_a + rho_b)
    end

    return outbound_dist
end

"""
BernoulliNode:

    Beta      Bern
    --->[Bern]--->
    <--
"""
function variationalRule!(  node::BernoulliNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Beta,
                            marg_in::Any,
                            marg_out::Bernoulli)

    outbound_dist.a = marg_out.p + 1.0
    outbound_dist.b = 2.0 - marg_out.p

    return outbound_dist
end


############################
# Average energy functional
############################

function U(::Type{BernoulliNode}, marg_in::Beta, marg_out::Bernoulli)
    digamma(marg_in.a + marg_in.b) -
    (1 - marg_out.p)*digamma(marg_in.b) -
    marg_out.p*digamma(marg_in.a)
end
