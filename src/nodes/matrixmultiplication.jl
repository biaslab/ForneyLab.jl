############################################
# MatrixMulitiplicationNode
############################################
# Description:
#   Multiplication with a predefined matrix:
#   out = A * in1
#   Example:
#       MatrixMulitiplicationNode([1.0]; name="my_node")
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

type MatrixMultiplicationNode <: Node
    A::Array
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    out::Interface
    function MatrixMultiplicationNode(A::Array; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        self = new(deepcopy(A), name, Array(Interface, 2))
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.out = self.interfaces[2]
        return self
    end
end
MatrixMultiplicationNode(; args...) = MatrixMultiplicationNode([1.0]; args...)

function calculatemessage!( outbound_interface_id::Int,
                            node::MatrixMultiplicationNode,
                            inbound_messages::Array{GaussianMessage,1})
    if outbound_interface_id == 1
        # Backward message
        msg_in = inbound_messages[2]
        msg = deepcopy(msg_in)
        if typeof(msg.xi)==Array
            msg.xi = transpose(node.A) * msg_in.xi
            msg.m = [] # Invalidate m, no need to pre-calculate
        end
        if typeof(msg.m)==Array || is(msg.W, nothing)
            if isposdef(node.A) && isposdef(msg_in.W)
                msg.m = inv(node.A) * msg_in.m
            else
                msg.m = pinv(transpose(node.A) * msg_in.W * node.A) * transpose(node.A) * msg_in.W * msg_in.m
            end
        end
        if typeof(msg.W)==Array
            msg.W = transpose(node.A) * msg_in.W * node.A
            msg.V = [] # Invalidate V, no need to pre-calculate
        end
        if typeof(msg.V)==Array
            # TODO: validate that A is non-singular
            A_inv = inv(node.A)
            msg.V = A_inv * msg_in.V * transpose(A_inv)
        end
    elseif outbound_interface_id == 2
        # Forward message
        msg_in = inbound_messages[1]
        msg = deepcopy(msg_in)
        if typeof(msg.m)==Array
            msg.m = node.A * msg_in.m
            msg.xi = nothing # Invalidate xi, no need to pre-calculate
        end
        if typeof(msg.xi)==Array
            if isposdef(node.A) && isposdef(msg_in.V)
                msg.xi = transpose(inv(node.A)) * msg_in.xi
            else
                msg.xi = pinv(node.A * msg_in.V * transpose(node.A)) * node.A * msg_in.V * msg_in.xi
            end
        end
        if typeof(msg.V)==Array
            msg.V = node.A * msg_in.V * transpose(node.A)
            msg.W = nothing # Invalidate W
        end
        if typeof(msg.W)==Array || typeof(msg.xi)==Array
            # Calculate W if V is unknown or xi is used
            # TODO: check if A is indeed non-singular and msg_in.W is available
            # TODO: otherwise, check if we can use W = inv(V) instead
            A_inv = inv(node.A)
            msg.W = transpose(A_inv) * msg_in.W * A_inv
        end
    end

    node.interfaces[outbound_interface_id].message = msg
end

function calculatemessage!( outbound_interface_id::Int,
                            node::MatrixMultiplicationNode,
                            inbound_messages::Array{GaussianMessage,1})
    if outbound_interface_id == 1
        # Backward message
        msg_in = inbound_messages[2].value
        msg = deepcopy(msg_in)
        if isposdef(node.A)
            msg.value = inv(node.A) * msg_in.value
        else
            msg.value = pinv(node.A) * msg_in.value
        end
    elseif outbound_interface_id == 2
        # Forward message
        msg = GeneralMessage(node.A * inbound_messages[1].value)
    end

    node.interfaces[outbound_interface_id].message = msg
end