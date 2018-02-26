export
ruleSPNonlinearOutVGP,
ruleSPNonlinearIn1GVP,
ruleVBNonlinearOut,
ruleVBNonlinearIn1,
ruleVBNonlinearW,
ruleSVBNonlinearOutVGD,
ruleSVBNonlinearW,
ruleSVBNonlinearIn1GVD,
ruleMNonlinearGGD


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

function ruleSPNonlinearIn1GVP{V<:VariateType}( msg_out::Message{Gaussian, V},
                                                msg_in1::Message, # Any type of message, used for determining approximation point
                                                msg_w::Message{PointMass},
                                                g::Function,
                                                J_g::Function)

    ensureParameters!(msg_out.dist, (:m, :v))
    (A, b) = approximate(unsafeMean(msg_in1.dist), g, J_g)
    A_inv = pinv(A)
    
    Message(V, Gaussian, m=A_inv*(msg_out.dist.params[:m] - b), v=A_inv*(msg_out.dist.params[:v] + cholinv(msg_w.dist.params[:m]))*A_inv')
end

function ruleVBNonlinearOut{V<:VariateType}(dist_out::Any,
                                            dist_in1::ProbabilityDistribution{V},
                                            dist_w::ProbabilityDistribution,
                                            g::Function,
                                            J_g::Function)

    (A, b) = approximate(unsafeMean(dist_in1), g, J_g)

    Message(V, Gaussian, m=A*unsafeMean(dist_in1) + b, w=unsafeMean(dist_w))
end

function ruleVBNonlinearIn1{V<:VariateType}(dist_out::ProbabilityDistribution{V},
                                            dist_in1::ProbabilityDistribution{V},
                                            dist_w::ProbabilityDistribution,
                                            g::Function,
                                            J_g::Function)

    (A, b) = approximate(unsafeMean(dist_in1), g, J_g) # Approximation is made on the marginal of in1

    Message(V, Gaussian, m=pinv(A)*(unsafeMean(dist_out) - b), w=A'*unsafeMean(dist_w)*A)
end

function ruleVBNonlinearW(  dist_out::ProbabilityDistribution{Univariate},
                            dist_in1::ProbabilityDistribution{Univariate},
                            dist_w::Any,
                            g::Function,
                            J_g::Function)

    (A, b) = approximate(unsafeMean(dist_in1), g, J_g)

    Message(Univariate, Gamma, a=1.5, b=0.5*(A^2*unsafeVar(dist_in1) + unsafeVar(dist_out) + (A*unsafeMean(dist_in1) + b - unsafeMean(dist_out))^2))
end

function ruleVBNonlinearW(  dist_out::ProbabilityDistribution{Multivariate},
                            dist_in1::ProbabilityDistribution{Multivariate},
                            dist_w::Any,
                            g::Function,
                            J_g::Function)

    (A, b) = approximate(unsafeMean(dist_in1), g, J_g)

    Message(MatrixVariate, Wishart, v=cholinv( A*unsafeCov(dist_in1)*A' + unsafeCov(dist_out) + (A*unsafeMean(dist_in1) + b - unsafeMean(dist_out))*(A*unsafeMean(dist_in1) + b - unsafeMean(dist_out))' ), nu=dims(dist_out) + 2.0)
end

function ruleSVBNonlinearOutVGD{V<:VariateType}(msg_out::Any,
                                                msg_in1::Message{Gaussian, V},
                                                dist_prec::ProbabilityDistribution,
                                                g::Function,
                                                J_g::Function)
    ensureParameters!(msg_in1.dist, (:m, :v))
    (A, b) = approximate(msg_in1.dist.params[:m], g, J_g)

    Message(V, Gaussian, m=A*msg_in1.dist.params[:m] + b, v=A*msg_in1.dist.params[:v]*A' + cholinv(unsafeMean(dist_prec)))
end

function ruleSVBNonlinearW( dist_out_in1::ProbabilityDistribution{Multivariate, Gaussian},
                            dist_prec::Any,
                            g::Function,
                            J_g::Function)

    ensureParameters!(dist_out_in1, (:m, :v))

    joint_dims = dims(dist_out_in1)
    if joint_dims == 2
        V = dist_out_in1.params[:v]
        m = dist_out_in1.params[:m]
        (A, b) = approximate(m[2], g, J_g)
        
        return Message(Univariate, Gamma, a=1.5, b=0.5*(A^2*V[2,2] - A*(V[1,2] + V[2,1]) + V[1,1] + (A*m[2] + b - m[1])^2))
    else
        V = dist_out_in1.params[:v]
        m = dist_out_in1.params[:m]
        d = Int64(joint_dims/2)
        (A, b) = approximate(m[d+1:end], g, J_g)
        
        return Message(MatrixVariate, Wishart, v=cholinv( A*V[d+1:end,d+1:end]*A' - A*V[d+1:end, 1:d] - V[1:d,d+1:end]*A' + V[1:d,1:d] + (A*m[d+1:end] + b - m[1:d])*(A*m[d+1:end] + b - m[1:d])' ), nu=d + 2.0) 
    end
end

function ruleSVBNonlinearIn1GVD{V<:VariateType}(msg_out::Message{Gaussian, V},
                                                msg_in1::Message, # Any type of message, used for determining approximation point
                                                dist_prec::ProbabilityDistribution,
                                                g::Function,
                                                J_g::Function)

    ensureParameters!(msg_out.dist, (:m, :v))
    (A, b) = approximate(unsafeMean(msg_in1.dist), g, J_g)
    A_inv = pinv(A)
    
    Message(V, Gaussian, m=A_inv*(msg_out.dist.params[:m] - b), v=A_inv*(msg_out.dist.params[:v] + cholinv(unsafeMean(dist_prec)))*A_inv')
end

function ruleMNonlinearGGD{V<:VariateType}( msg_out::Message{Gaussian, V},
                                            msg_in1::Message{Gaussian, V},
                                            dist_prec::ProbabilityDistribution,
                                            g::Function,
                                            J_g::Function)

    ensureParameters!(msg_out.dist, (:m, :w))
    ensureParameters!(msg_in1.dist, (:m, :w))
    (A, b) = approximate(msg_in1.dist.params[:m], g, J_g)

    m_y = msg_out.dist.params[:m]
    m_x = msg_in1.dist.params[:m]
    W_y = msg_out.dist.params[:w]
    W_x = msg_in1.dist.params[:w]
    W_bar = unsafeMean(dist_prec)

    V_q = cholinv([W_y+W_bar -W_bar*A; -A'*W_bar W_x+A'*W_bar*A])

    return ProbabilityDistribution(Multivariate, Gaussian, m=V_q*[W_y*m_y+W_bar*b; W_x*m_x-A'*W_bar*b], v=V_q)
end


#--------------------------
# Custom inbound collectors
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

function collectNaiveVariationalNodeInbounds(node::Nonlinear, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        if node_interface == node.interfaces[2]
            # Always collect the marginal on the in1 edge,
            # because it is used to obtain the approximation point
            push!(inbounds, "marginals[:$(node_interface.edge.variable.id)]")
        elseif node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        else
            # Collect marginal from marginal dictionary
            push!(inbounds, "marginals[:$(node_interface.edge.variable.id)]")
        end
    end

    # Push function and Jacobi matrix function names to calling signature
    # These functions need to be defined in the scope of the user
    push!(inbound_messages, "$(node.g)")
    push!(inbound_messages, "$(node.J_g)")

    return inbounds
end

function collectStructuredVariationalNodeInbounds(node::Nonlinear, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(entry.interface.edge)
    local_cluster_ids = localRecognitionFactorization(entry.interface.node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if (node_interface == node.interfaces[2]) && (node_interface_recognition_factor_id == entry_recognition_factor_id)
            # Always collect the message inbound on the in1 edge,
            # because it is used to obtain the approximation point
            haskey(interface_to_msg_idx, inbound_interface) || error("The nonlinear node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        elseif node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    # Push function and Jacobi matrix function names to calling signature
    # These functions need to be defined in the scope of the user
    push!(inbounds, "$(node.g)")
    push!(inbounds, "$(node.J_g)")

    return inbounds
end

function collectMarginalNodeInbounds(node::Nonlinear, entry::MarginalScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(first(entry.target.edges))
    local_cluster_ids = localRecognitionFactorization(node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    # Push function and Jacobi matrix function names to calling signature
    # These functions need to be defined in the scope of the user
    push!(inbounds, "$(node.g)")
    push!(inbounds, "$(node.J_g)")

    return inbounds
end

function collectAverageEnergyInbounds(node::Nonlinear)
    inbounds = String[]

    local_cluster_ids = localRecognitionFactorization(node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    # Push function and Jacobi matrix function names to calling signature
    # These functions need to be defined in the scope of the user
    push!(inbounds, "$(node.g)")
    push!(inbounds, "$(node.J_g)")

    return inbounds
end
