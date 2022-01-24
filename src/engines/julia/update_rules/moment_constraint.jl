export ruleSPMomentConstraintOutG

"""
Evaluate the expectation ∫q(x; η)g(x)dx for a given eta
"""
function expectedValue(eta::Float64, m_bw::Float64, V_bw::Float64, g::Function)
    Z_tilde = gaussianQuadrature(x->exp(eta*g(x)),              m=m_bw, v=V_bw) # Obtain normalizing constant for constrained q(x)
    G_tilde = gaussianQuadrature(x->exp(eta*g(x))*g(x)/Z_tilde, m=m_bw, v=V_bw) # Evaluate the expectation of g(x) w.r.t. the constrained q(x)

    return G_tilde
end

"""
Obtain the first and second moment of the constrained q(x; η) for a given η
"""
function constrainedMoments(eta::Float64, m_bw::Float64, V_bw::Float64, g::Function)
    Z_tilde    = gaussianQuadrature(x->exp(eta*g(x)),             m=m_bw, v=V_bw) # Obtain normalizing constant for constrained q(x)
    m_tilde    = gaussianQuadrature(x->exp(eta*g(x))*x/Z_tilde,   m=m_bw, v=V_bw) # Compute mean
    V_tilde_nc = gaussianQuadrature(x->exp(eta*g(x))*x^2/Z_tilde, m=m_bw, v=V_bw) # Compute non-central variance
    V_tilde    = V_tilde_nc - m_tilde^2 # Compute central variance

    return (m_tilde, V_tilde)
end


#-------------
# Update rules
#-------------

function ruleSPMomentConstraintOutG(msg_out::Message{<:Gaussian, Univariate}, g::Function, G::Float64, eta_init::Float64)
    (m_bw, V_bw) = unsafeMeanCov(msg_out.dist)
    (xi_bw, W_bw) = unsafeWeightedMeanPrecision(msg_out.dist)

    eta_star = find_zero(eta->expectedValue(eta, m_bw, V_bw, g)-G, eta_init) # Find optimal eta
    (m_tilde, V_tilde) = constrainedMoments(eta_star, m_bw, V_bw, g) # Compute moments of constrained q(x) with optimal eta
    
    # Convert approximated marginal to canonical form
    W_tilde = inv(V_tilde)
    xi_tilde = W_tilde*m_tilde

    return Message(Univariate, Gaussian{Canonical}, xi=xi_tilde-xi_bw, w=W_tilde-W_bw) # Compute forward message
end


#---------------------------
# Custom inbounds collectors
#---------------------------

function collectSumProductNodeInbounds(node::MomentConstraint, entry::ScheduleEntry)
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
    push!(inbounds, Dict{Symbol, Any}(:G => node.G,
                                      :keyword => false))
    push!(inbounds, Dict{Symbol, Any}(:eta_init => node.eta_init,
                                      :keyword => false))

    return inbounds
end