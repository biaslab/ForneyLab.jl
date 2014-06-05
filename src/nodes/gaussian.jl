############################################
# GaussianNode
############################################
# Description:
#   Node converting an input mean and precision
#   to a univariate Gaussian distribution:
#
#         in1 (mean)
#          |
#          v  out
#   ----->[N]----->
#  in2 (prec.)
#
#   out = GaussianMessage(m=in1, W=in2)
#
#   Example:
#       GaussianNode([1.0], [0.1]; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       GaussianMessage
#   2. (in2):
#       GammaMessage
#   3. (out):
#       GaussianMessage
############################################

export GaussianNode

type GaussianNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface
    function GaussianNode(; args...)
        (name = getArgumentValue(args, :name))!=false || (name = "unnamed")
        self = new(name, Array(Interface, 3))
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.in2 = self.interfaces[2]
        self.out = self.interfaces[3]
        return self
    end
end

############################################
# Gaussian- and GammaMessage methods
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::GaussianNode,
                            inbound_messages_types::Type{Union(GaussianMessage, GammaMessage)})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    if outbound_interface_id == 3
        # Forward message
        mean_msg = node.interfaces[1].partner.message
        prec_msg = node.interfaces[2].partner.message
        if typeof(mean_msg) != GaussianMessage || typeof(prec_msg) != GaussianMessage
            error("Incompatible message type with update function for node $(typeof(node)) $(node.name).")
        end
        # TODO: now we take the expectations of the incoming messages, do the variational approach
        msg_out = GaussianMessage(m = mean_msg.m, W=prec_msg.a/prec_msg.b)
    else
        # Backward message
        # TODO: implement
        error("Backward messages not implemented yet for $(typeof(node)) type.")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg_out
end

