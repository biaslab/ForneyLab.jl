export
ruleSPNonlinearOutVGP,
ruleVBNonlinearOut,
ruleVBNonlinearIn1,
ruleVBNonlinearW

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
    push!(inbound_messages, "$(node.g)")
    push!(inbound_messages, "$(node.J_g)")

    return inbounds
end
