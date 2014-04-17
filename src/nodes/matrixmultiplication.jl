############################################
# MatrixMulitiplicationNode
############################################
# Description:
#   Multiplication with a predefined value or matrix
#   out = A * in1
#   One can define the type of the matrix values
#   i.e. MatrixMulitiplicationNode{Float64}("unity_node", 1.0)
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       GaussianMessage
#       GeneralMessage
#   2. (out):
#       GaussianMessage
#       GeneralMessage
############################################

export MatrixMultiplicationNode

type MatrixMultiplicationNode{T} <: Node
    interfaces::Array{Interface,1}
    name::ASCIIString
    A::Array{T,2}
    in1::Interface
    out::Interface
    function MatrixMultiplicationNode(name::ASCIIString, A::Array)
        self = new(Array(Interface, 2), name)
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.out = self.interfaces[3]
        # Store A
        self.A = deepcopy(A)
        return self
    end
    function MatrixMultiplicationNode(A::Array)
        return MatrixMultiplicationNode("#undef", A)
    end
end

function calculatemessage!{T<:Union(GaussianMessage,GeneralMessage)}(
                            outbound_interface_id::Int,
                            node::MatrixMultiplicationNode,
                            inbound_messages::Array{T,1})
    if outbound_interface_id == 1 # message to in1 interface (backward multiplicaton)
        # TODO
        msg_in = inbound_messages[2]
        msg = deepcopy(msg_in)
        node.interfaces[outbound_interface_id].message = msg
    elseif outbound_interface_id == 2 # message to out interface (forward multiplicaton)
        if typeof(inbound_messages[1]) == GeneralMessage
            node.interfaces[outbound_interface_id].message = GeneralMessage(node.A * inbound_messages[1].value)
        elseif typeof(inbound_messages[1]) == GaussianMessage
            msg_in = inbound_messages[1]
            msg = deepcopy(msg_in)
            if length(msg.m)>0
                msg.m = node.A * msg_in.m
                msg.xi = [] # Invalidate xi, no need to precalculate it
            end
            if length(msg.xi)>0
                if isposdef(node.A) && isposdef(msg_in.V)
                    msg.xi = transpose(node.A) * msg_in.xi
                else
                    msg.xi = pinv(node.A * msg_in.V * transpose(node.A)) * node.A * msg_in.V * msg_in.xi
                end
            end
            if length(msg.V)>0
                msg.V = node.A * msg_in.V * transpose(node.A)
                msg.W = [] # Invalidate W, no need to precalculate it
            end
            if length(msg.W)>0
                A_inv = inv(node.A) # NOTE: we assume that A is non-singular here
                # TODO: see if it is beneficial to cache A_inv (compiler might already do this)
                msg.W = transpose(A_inv) * msg_in.V * A_inv
            end
            node.interfaces[outbound_interface_id].message = msg
        end
    end
end