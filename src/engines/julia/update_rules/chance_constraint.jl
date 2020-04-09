export ruleSPChanceConstraintOutG

# TODO: implement multivariate update using cubature
# TODO: enforce g(x) ∈ {0, 1}
function ruleSPChanceConstraintOutG(msg_out::Message{<:Gaussian, Univariate}, g::Function, epsilon::Float64)
    (m_bw, V_bw) = unsafeMeanCov(msg_out.dist)
    (xi_bw, W_bw) = unsafeWeightedMeanPrecision(msg_out.dist)

    # Default forward message is vague
    xi_fw = 0.0
    W_fw = tiny

    # Correct marginal belief (if required)
    Phi_G = gaussianQuadrature(g, m=m_bw, v=V_bw)
    if epsilon <= 1.0 - Phi_G # If constraint is active
        # Non-central moments of individual regions
        m_G     = 1.0/Phi_G*gaussianQuadrature(x->g(x)*x, m=m_bw, v=V_bw)
        m_nG    = 1.0/(1.0 - Phi_G)*gaussianQuadrature(x->(1.0 - g(x))*x, m=m_bw, v=V_bw)
        V_nc_G  = 1.0/Phi_G*gaussianQuadrature(x->g(x)*x^2, m=m_bw, v=V_bw)
        V_nc_nG = 1.0/(1.0 - Phi_G)*gaussianQuadrature(x->(1.0 - g(x))*x^2, m=m_bw, v=V_bw)

        # Central moments of corrected belief
        m_tilde = (1 - epsilon)*m_G + epsilon*m_nG # Mean
        V_nc_tilde = (1 - epsilon)*V_nc_G + epsilon*V_nc_nG # Noncentral variance
        V_tilde = V_nc_tilde - m_tilde^2 # Central variance

        # Convert approximated marginal to canonical form
        W_tilde = inv(V_tilde)
        xi_tilde = W_tilde*m_tilde

        # Compute forward message
        xi_fw = xi_tilde - xi_bw
        W_fw  = W_tilde - W_bw
    end

    return Message(Univariate, GaussianWeightedMeanPrecision, xi=xi_fw, w=W_fw)
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::ChanceConstraint, entry::ScheduleEntry)
    inbounds = Any[]

    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    # Push function and value to calling signature; these need to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))
    push!(inbounds, Dict{Symbol, Any}(:epsilon => node.epsilon,
                                      :keyword => false))

    return inbounds
end