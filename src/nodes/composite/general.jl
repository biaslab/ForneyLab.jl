# Implementation of updateNodeMessage!() for a general CompositeNode with an internal graph (use_composite_update_rules==false).
# This implementation can be used if the outbound messages can be calculated by internal message passing.
# The internal graph should have no loops. If it does, the user should overload updateNodeMessage!() to set initial messages, iterate through the loop(s) and check convergence.

function updateNodeMessage!(outbound_interface_id::Int,
                            node::CompositeNode,
                            inbound_messages_types::Type{GaussianMessage})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    if node.use_composite_update_rules
    	error("$(typeof(node)) $(node.name) is configured to use shortcut rules, but updateNodeMessage!() is not defined for this node type.")
	else
		# Internal message passing
	    msg_out = calculateMessage!(node.interfaces[outbound_interface_id].child)
	    # Copy the outbound message to composite node's interface and return it
	    return node.interfaces[outbound_interface_id].message = msg_out
	end
end