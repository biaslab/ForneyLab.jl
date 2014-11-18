# Implementation of updateNodeMessage!() for a general CompositeNode with an internal graph (use_composite_update_rules==false).
# This implementation can be used if the outbound messages can be calculated by internal message passing.
# The internal graph should have no loops. If it does, the user should overload updateNodeMessage!() to set initial messages, iterate through the loop(s) and check convergence.

function updateNodeMessage!(node::CompositeNode,
                            outbound_interface_id::Int,
                            inbounds...)
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    if node.use_composite_update_rules
        error("$(typeof(node)) $(node.name) is configured to use shortcut rules, but updateNodeMessage!() with arguments $([typeof(inbound) for inbound in [inbounds...]]) is not defined for this node type.")
    else
        # Internal message passing
        schedule = nothing
        try
            schedule = node.interfaces[outbound_interface_id].internal_schedule
        catch
            error("No internal message passing schedule is defined for calculating outbound messages on interface $(outbound_interface_id) of $(typeof(node)) $(node.name).")
        end

        return node.interfaces[outbound_interface_id].message = executeSchedule(schedule)

    end
end