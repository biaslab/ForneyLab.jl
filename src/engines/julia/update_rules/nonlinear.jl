export
ruleSPNonlinearOutVG,
ruleSPNonlinearIn1GV


"""
Find a local linear approximation to the nonlinear vector function g at x_hat
"""
function approximate(x_hat::Union{Float64, Vector{Float64}}, g::Function, J_g::Function)
    A = J_g(x_hat)
    b = g(x_hat) - A*x_hat

    return (A, b)
end

function ruleSPNonlinearOutVG{V<:VariateType}(  msg_out::Void,
                                                msg_in1::Message{Gaussian, V},
                                                g::Function,
                                                J_g::Function)

    ensureParameters!(msg_in1.dist, (:m, :v))
    (A, b) = approximate(msg_in1.dist.params[:m], g, J_g)
    (A != 0.0) || (A = tiny*sign(randn())) # Perturb A

    Message(V, Gaussian, m=A*msg_in1.dist.params[:m] + b, v=A*msg_in1.dist.params[:v]*A')
end

function ruleSPNonlinearIn1GV{V<:VariateType}(  msg_out::Message{Gaussian, V},
                                                msg_in1::Message, # Any type of message, used for determining the approximation point
                                                g::Function,
                                                J_g::Function)

    ensureParameters!(msg_out.dist, (:m, :w))
    (A, b) = approximate(unsafeMean(msg_in1.dist), g, J_g)
    (A != 0.0) || (A = tiny*sign(randn())) # Perturb A
    A_inv = pinv(A)
    
    Message(V, Gaussian, m=A_inv*(msg_out.dist.params[:m] - b), w=A'*msg_out.dist.params[:w]*A)
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::Nonlinear, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == node.interfaces[2]
            # Always collect the message inbound on the in1 edge,
            # because it is used to obtain the approximation point
            haskey(interface_to_msg_idx, inbound_interface) || error("The nonlinear node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        elseif node_interface == entry.interface
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
