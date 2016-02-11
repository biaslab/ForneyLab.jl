export TerminalNode, PriorNode

"""
Description:

    Sends out a predefined message.

        out
    [T]----->

    out = T.value

Interfaces:
    
    1 i[:out]

Construction:
    
    TerminalNode(GaussianDistribution(), id=:my_node)
"""
type TerminalNode <: Node
    id::Symbol
    value::ProbabilityDistribution
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function TerminalNode(value::ProbabilityDistribution; id=generateNodeId(TerminalNode))
        self = new(id, deepcopy(value), Vector{Interface}(1), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        self.i[:out] = self.interfaces[1] = Interface(self)

        return self
    end
end

TerminalNode(num::Number; id=generateNodeId(TerminalNode)) = TerminalNode(convert(DeltaDistribution, num), id=id)

TerminalNode(bool::Bool; id=generateNodeId(TerminalNode)) = TerminalNode(convert(DeltaDistribution, bool), id=id)

TerminalNode{T<:Number}(vect::Vector{T}; id=generateNodeId(TerminalNode)) = TerminalNode(convert(MvDeltaDistribution, vect), id=id)

TerminalNode(; id=generateNodeId(TerminalNode)) = TerminalNode(vague(GaussianDistribution), id=id)

typealias PriorNode TerminalNode # For more overview during graph construction

isDeterministic(::TerminalNode) = false # Edge case for deterministicness

# Implement firstFreeInterface since EqualityNode is symmetrical in its interfaces
firstFreeInterface(node::TerminalNode) = (node.interfaces[1].partner==nothing) ? node.interfaces[1] : error("No free interface on $(typeof(node)) $(node.id)")

function sumProductRule!(   node::TerminalNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Any,
                            msg_out::Any)

    # Fill the fields of outbound_dist with node.value
    return injectParameters!(outbound_dist, node.value)
end
