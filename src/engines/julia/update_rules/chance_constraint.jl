export ruleSPChanceConstraintOutG

const default_atol = 1e-4

"""
Computes normalizing constant and central moments of a truncated Gaussian distribution.
"""
function truncatedGaussianMoments(m::Float64, V::Float64, a::Float64, b::Float64)
    # Normalize parameters
    sigma = sqrt(clamp(V, tiny, huge))
    alpha = (a - m)/sigma
    beta = (b - m)/sigma
    
    # Compute normalizing constant
    Phi_alpha = normcdf(alpha)
    Phi_beta = normcdf(beta)
    Z = Phi_beta - Phi_alpha

    # Compute moments if truncation is valid
    if Z < tiny
        # Invalid region; return undefined mean and variance of truncated distribution
        Z    = 0.0
        m_tr = 0.0
        V_tr = 0.0
    else
        phi_alpha = normpdf(alpha)
        phi_beta = normpdf(beta)
        alpha_phi_alpha = clamp(alpha, -huge, huge)*phi_alpha
        beta_phi_beta = clamp(beta, -huge, huge)*phi_beta
    
        m_tr = m + (phi_alpha - phi_beta)/Z*sigma
        V_tr = V*( 1.0 + (alpha_phi_alpha - beta_phi_beta)/Z - ((phi_alpha - phi_beta)/Z)^2 )
    end

    return (Z, m_tr, V_tr)
end


#-------------
# Update rules
#-------------

function ruleSPChanceConstraintOutG(msg_out::Message{<:Gaussian, Univariate}, G::Tuple, epsilon::Float64; atol=default_atol)
    (m_bw, V_bw) = unsafeMeanCov(msg_out.dist)
    (xi_bw, W_bw) = unsafeWeightedMeanPrecision(msg_out.dist)

    # Default forward message is vague; only corrected when constraint is active
    xi_fw = 0.0
    W_fw = tiny

    # Correct marginal belief if constraint is active (1-D only)
    (min_G, max_G) = G # Extract minimum and maximum of safe region G
    (Phi_G, m_G, V_G) = truncatedGaussianMoments(m_bw, V_bw, min_G, max_G) # Compute statistics (and normalizing constant) of q in G
    if epsilon <= 1.0 - Phi_G # If constraint is active
        # Initialize statistics of uncorrected belief
        m_tilde = m_bw
        V_tilde = V_bw
        for i = 1:100 # Iterate at most this many times
            (Phi_lG, m_lG, V_lG) = truncatedGaussianMoments(m_tilde, V_tilde, -Inf, min_G) # Statistics for q in region left of G
            (Phi_rG, m_rG, V_rG) = truncatedGaussianMoments(m_tilde, V_tilde, max_G, Inf) # Statistics for q in region right of G

            # Compute moments of non-G region as a mixture of left and right truncations
            Phi_nG = Phi_lG + Phi_rG
            m_nG = Phi_lG/Phi_nG*m_lG + Phi_rG/Phi_nG*m_rG
            V_nG = Phi_lG/Phi_nG*(V_lG + m_lG^2) + Phi_rG/Phi_nG*(V_rG + m_rG^2) - m_nG^2

            # Compute moments of corrected belief as a mixture of G and non-G regions
            m_tilde = (1.0-epsilon)*m_G + epsilon*m_nG
            V_tilde = (1.0-epsilon)*(V_G + m_G^2) + epsilon*(V_nG + m_nG^2) - m_tilde^2

            # Re-compute statistics (and normalizing constant) of corrected belief
            (Phi_G, m_G, V_G) = truncatedGaussianMoments(m_tilde, V_tilde, min_G, max_G)

            if (1.0 - Phi_G) < (1.0 + atol)*epsilon
                break # Break the loop if the belief is sufficiently corrected
            end
        end

        # Convert moments of corrected belief to canonical form
        W_tilde = cholinv(V_tilde)
        xi_tilde = W_tilde*m_tilde

        # Compute canonical parameters of forward message
        xi_fw = xi_tilde - xi_bw
        W_fw  = W_tilde - W_bw
    end

    return Message(Univariate, Gaussian{Canonical}, xi=xi_fw, w=W_fw)
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
    push!(inbounds, Dict{Symbol, Any}(:G => node.G,
                                      :keyword => false))
    push!(inbounds, Dict{Symbol, Any}(:epsilon => node.epsilon,
                                      :keyword => false))
    if node.atol != nothing
        push!(inbounds, Dict{Symbol, Any}(:atol => node.atol,
                                          :keyword => true))
    end

    return inbounds
end