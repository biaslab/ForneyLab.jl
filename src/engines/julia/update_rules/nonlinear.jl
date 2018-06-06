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

function ruleSPNonlinearOutVG(  msg_out::Void,
                                msg_in1::Message{Gaussian, Multivariate},
                                g::Function,
                                J_g::Function)

    ensureParameters!(msg_in1.dist, (:m, :v))
    (A, b) = approximate(msg_in1.dist.params[:m], g, J_g)
    V_q = A*msg_in1.dist.params[:v]*A'
    V_q = V_q + tiny*diageye(size(V_q)[1]) # Ensure V_q is invertible

    Message(Multivariate, Gaussian, m=A*msg_in1.dist.params[:m] + b, v=V_q)
end

function ruleSPNonlinearOutVG(  msg_out::Void,
                                msg_in1::Message{Gaussian, Univariate},
                                g::Function,
                                J_g::Function)

    ensureParameters!(msg_in1.dist, (:m, :v))
    (a, b) = approximate(msg_in1.dist.params[:m], g, J_g)
    v_q = clamp(msg_in1.dist.params[:v]*a^2, tiny, huge)

    Message(Univariate, Gaussian, m=a*msg_in1.dist.params[:m] + b, v=v_q)
end

function ruleSPNonlinearIn1GV(  msg_out::Message{Gaussian, Multivariate},
                                msg_in1::Message, # Any type of message, used for determining the approximation point
                                g::Function,
                                J_g::Function)

    ensureParameters!(msg_out.dist, (:m, :w))
    (A, b) = approximate(unsafeMean(msg_in1.dist), g, J_g)
    A_inv = pinv(A)
    W_q = A'*msg_out.dist.params[:w]*A
    W_q = W_q + tiny*diageye(size(W_q)[1]) # Ensure W_q is invertible
    
    Message(Multivariate, Gaussian, m=A_inv*(msg_out.dist.params[:m] - b), w=W_q)
end

function ruleSPNonlinearIn1GV(  msg_out::Message{Gaussian, Univariate},
                                msg_in1::Message, # Any type of message, used for determining the approximation point
                                g::Function,
                                J_g::Function)

    ensureParameters!(msg_out.dist, (:m, :w))
    (a, b) = approximate(unsafeMean(msg_in1.dist), g, J_g)
    w_q = clamp(msg_out.dist.params[:w]*a^2, tiny, huge)
    
    Message(Univariate, Gaussian, m=(msg_out.dist.params[:m] - b)/a, w=w_q)
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
