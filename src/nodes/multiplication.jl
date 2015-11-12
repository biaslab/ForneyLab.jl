############################################
# MultiplicationNode
############################################
# Description:
#   Multiplication: out = in1 * in2
#
#          in2
#          |
#    in1   v  out
#   ----->[+]----->
#
#   f(in1,in2,out) = δ(out - (in1 * in2))
#
# Interfaces:
#   1 i[:in1], 2 i[:in2], 3 i[:out]
#
# Construction:
#   MultiplicationNode(id=:my_node)
#
############################################

export MultiplicationNode

type MultiplicationNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function AdditionNode(; id=generateNodeId(MultiplicationNode))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        addNode!(current_graph, self)

        for (iface_index, iface_handle) in enumerate([:in1, :in2, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::MultiplicationNode) = true

###################################################
### Gaussian and Delta methods
###################################################
