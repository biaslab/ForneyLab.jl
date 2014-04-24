############################################
# EqualityNode
############################################
# Description:
#   Equality constraint node with variable
#   number of (symmetrical) interfaces (ports).
#   Example:
#       EqualityNode(name="3_port_equ") # 3 interfaces is default
#       EqualityNode(5; name="5_port_equ")
# Interface ids, (names) and supported message types:
#   1. (none):
#       GaussianMessage
#   2. (none):
#       GaussianMessage
#   3. (none):
#       GaussianMessage
#   ...
#   N. (none):
#       GaussianMessage
############################################

export EqualityNode

type EqualityNode <: Node
    num_interfaces::Uint16
    interfaces::Array{Interface,1}
    name::ASCIIString

    function EqualityNode(num_interfaces::Integer; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        @assert(num_interfaces>2, "An EqualityNode should have at least 3 interfaces")
        self = new(num_interfaces, Array(Interface, num_interfaces), name)
        # Create interfaces
        for interface_id=1:num_interfaces
            self.interfaces[interface_id] = Interface(self)
        end
        return self
    end
end
EqualityNode(; args...) = EqualityNode(3; args...)

############################################
# GaussianMessage methods
############################################

function calculateMessage!( outbound_interface_id::Int,
                            node::EqualityNode,
                            inbound_messages::Array{GaussianMessage,1})
    # TODO: implement
    return node.interfaces[outbound_interface_id].message = GaussianMessage()
end