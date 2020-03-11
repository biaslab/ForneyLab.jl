export
ruleSPNonlinearGaussianOutNGP,
ruleSPNonlinearGaussianIn1GGP,
ruleMNonlinearGaussianGGD

# Forward rule
function ruleSPNonlinearGaussianOutNGP(msg_out::Nothing,
                                       msg_in1::Message{F, V},
                                       msg_var::Message{PointMass},
                                       g::Function;
                                       alpha::Float64=default_alpha) where {F<:Gaussian, V<:VariateType}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    (m_tilde, V_tilde, _) = unscentedStatistics(m_fw_in1, V_fw_in1, g; alpha=alpha)

    return Message(V, GaussianMeanVariance, m=m_tilde, v=V_tilde+msg_var.dist.params[:m]) # Include additive noise
end

# Backward rule with given inverse
function ruleSPNonlinearGaussianIn1GGP(msg_out::Message{F, V},
                                       msg_in1::Nothing,
                                       msg_var::Message{PointMass},
                                       g::Function,
                                       g_inv::Function;
                                       alpha::Float64=default_alpha) where {F<:Gaussian, V<:VariateType}
    
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    (m_tilde, V_tilde, _) = unscentedStatistics(m_bw_out, V_bw_out+msg_var.dist.params[:m], g_inv; alpha=alpha) # Include additive noise

    return Message(V, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

# Backward rule with unknown inverse
function ruleSPNonlinearGaussianIn1GGP(msg_out::Message{F1, V},
                                       msg_in1::Message{F2, V},
                                       msg_var::Message{PointMass},
                                       g::Function;
                                       alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian, V<:VariateType}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    W_fw_in1 = unsafePrecision(msg_in1.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)

    (m_tilde, V_tilde, C_tilde) = unscentedStatistics(m_fw_in1, V_fw_in1, g; alpha=alpha)
    V_tilde += msg_var.dist.params[:m] # Include additive noise

    # Update based on (Petersen et al. 2018; On Approximate Nonlinear Gaussian Message Passing on Factor Graphs)
    # Note, this implementation is not as efficient as Petersen et al. (2018), because we explicitly require the outbound messages
    C_tilde_inv = pinv(C_tilde)
    V_bw_in1 = V_fw_in1*C_tilde_inv'*(V_tilde + V_bw_out)*C_tilde_inv*V_fw_in1 - V_fw_in1
    m_bw_in1 = m_fw_in1 - (V_fw_in1 + V_bw_in1)*W_fw_in1*C_tilde*cholinv(V_tilde + V_bw_out)*(m_tilde - m_bw_out)

    return Message(V, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end

# Marginal rule
function ruleMNonlinearGaussianGGD(msg_out::Message{F1, V},
                                   msg_in1::Message{F2, V},
                                   marg_var::ProbabilityDistribution,
                                   g::Function;
                                   alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian, V<:VariateType}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    (xi_fw_in1, W_fw_in1) = unsafeWeightedMeanPrecision(msg_in1.dist)
    (xi_bw_out, W_bw_out) = unsafeWeightedMeanPrecision(msg_out.dist)

    (m_tilde, V_tilde, C_tilde) = unscentedStatistics(m_fw_in1, V_fw_in1, g; alpha=alpha)
    V_tilde += unsafeMean(marg_var) # Include additive noise
    
    # Compute joint statistics
    z = m_tilde - C_tilde'*xi_fw_in1
    G_tilde = cholinv(V_tilde - C_tilde'*W_fw_in1*C_tilde)
    K = W_fw_in1*C_tilde*G_tilde
    xi_q_out = G_tilde*z + xi_bw_out
    xi_q_in1 = xi_fw_in1 - K*z
    W_q_out = G_tilde + W_bw_out
    W_q_in1 = W_fw_in1 + K*C_tilde'*W_fw_in1

    return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[xi_q_out; xi_q_in1], w=[W_q_out -K'; -K W_q_in1])
end