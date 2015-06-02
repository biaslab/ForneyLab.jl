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

type TerminalNode <: Node
    id::Symbol
    value::ProbabilityDistribution
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function TerminalNode(value=DeltaDistribution(1.0); id=generateNodeId(TerminalNode))
        if typeof(value) <: Message || typeof(value) == DataType
            error("TerminalNode $(id) can not hold value of type $(typeof(value)).")
        end
        self = new(id, deepcopy(value), Array(Interface, 1), Dict{Symbol,Interface}())
        !haskey(current_graph.n, id) ? current_graph.n[id] = self : error("Node id $(id) already present")

        self.i[:out] = self.interfaces[1] = Interface(self)
 
        return self
    end
end

typealias PriorNode TerminalNode # For more overview during graph construction

isDeterministic(::TerminalNode) = false # Edge case for deterministicness

# Implement firstFreeInterface since EqualityNode is symmetrical in its interfaces
firstFreeInterface(node::TerminalNode) = (node.interfaces[1].partner==nothing) ? node.interfaces[1] : error("No free interface on $(typeof(node)) $(node.id)")

function sumProduct!(node::TerminalNode,
                            outbound_interface_index::Int,
                            ::Any)
    # Calculate an outbound message. The TerminalNode does not accept incoming messages.
    # This function is not exported, and is only meant for internal use.
    if (typeof(node.interfaces[1].message) != Message{typeof(node.value)}) || (node.interfaces[1].message.payload != node.value)
        # Only create a new message if the existing one is not correct
        node.interfaces[1].message = Message(node.value)
    end

    return (:terminal_forward,
            node.interfaces[1].message)
end