export
ruleSPNonlinearOutVGP

"""
Find a local linear approximation to the nonlinear vector function g at x_hat
"""
function approximate(x_hat::Union{Float64, Vector{Float64}}, g::Function, J_g::Function)
    A = J_g(x_hat)
    b = g(x_hat) - A*x_hat

    return (A, b)
end

function ruleSPNonlinearOutVGP{V<:VariateType}( msg_out::Void,
                                                msg_in1::Message{Gaussian, V},
                                                msg_w::Message{PointMass},
                                                g::Function,
                                                J_g::Function)

    ensureParameters!(msg_in1.dist, (:m, :v))
    (A, b) = approximate(msg_in1.dist.params[:m], g, J_g)

    Message(V, Gaussian, m=A*msg_in1.dist.params[:m] + b, v=A*msg_in1.dist.params[:v]*A' + cholinv(msg_w.dist.params[:m]))
end

function collectSumProductNodeInbounds(node::Nonlinear, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbound_messages, "nothing")
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbound_messages, messageString(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        end
    end

    # Push function and Jacobi matrix function names to calling signature
    # These functions need to be defined in the scope of the user
    push!(inbound_messages, "$(node.g)")
    push!(inbound_messages, "$(node.J_g)")

    return inbound_messages
end