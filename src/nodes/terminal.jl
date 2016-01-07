############################################
# TerminalNode
############################################
# Description:
#   Sends out a predefined message.
#
#       out
#   [T]----->
#
#   out = T.value
#
# Interfaces:
#   1 i[:out]
#
# Construction:
#   TerminalNode(GaussianDistribution(), id=:my_node)
#
############################################

export TerminalNode, PriorNode

type TerminalNode{value_type<:ProbabilityDistribution} <: Node
    id::Symbol
    value::value_type
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function TerminalNode{T}(value::T=DeltaDistribution(); id=generateNodeId(TerminalNode))
        self = new(id, deepcopy(value), Vector{Interface}(1), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        self.i[:out] = self.interfaces[1] = Interface(self)

        return self
    end
end

TerminalNode(value::ProbabilityDistribution; id=generateNodeId(TerminalNode)) = TerminalNode{typeof(value)}(value, id=id)

TerminalNode(; id=generateNodeId(TerminalNode)) = TerminalNode{GaussianDistribution}(vague(GaussianDistribution), id=id)

typealias PriorNode TerminalNode # For more overview during graph construction

isDeterministic(::TerminalNode) = false # Edge case for deterministicness

# Implement firstFreeInterface since EqualityNode is symmetrical in its interfaces
firstFreeInterface(node::TerminalNode) = (node.interfaces[1].partner==nothing) ? node.interfaces[1] : error("No free interface on $(typeof(node)) $(node.id)")

function sumProduct!(   node::TerminalNode,
                        outbound_interface_index::Type{Val{1}},
                        msg_out::Any,
                        outbound_dist::Any)

    # Fill the fields of outbound_dist with node.value
    for field in fieldnames(outbound_dist)
        setfield!(outbound_dist, field, deepcopy(getfield(node.value, field)))
    end
    
    return outbound_dist
end