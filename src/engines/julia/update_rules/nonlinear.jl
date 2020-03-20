export
ruleSPNonlinearOutNG,
ruleSPNonlinearIn1GG,
ruleMNonlinearGG

const default_alpha = 1e-3 # Default value for the spread parameter
const default_beta = 2.0
const default_kappa = 0.0

"""
Return the sigma points and weights for a Gaussian distribution
"""
function sigmaPointsAndWeights(m::Float64, V::Float64; alpha=default_alpha, beta=default_beta, kappa=default_kappa)
    lambda = (1 + kappa)*alpha^2 - 1

    sigma_points = Vector{Float64}(undef, 3)
    weights_m = Vector{Float64}(undef, 3)
    weights_c = Vector{Float64}(undef, 3)

    l = sqrt((1 + lambda)*V)

    sigma_points[1] = m
    sigma_points[2] = m + l
    sigma_points[3] = m - l
    weights_m[1] = lambda/(1 + lambda)
    weights_m[2] = weights_m[3] = 1/(2*(1 + lambda))
    weights_c[1] = weights_m[1] + (1 - alpha^2 + beta)
    weights_c[2] = weights_c[3] = 1/(2*(1 + lambda))

    return (sigma_points, weights_m, weights_c)
end

function sigmaPointsAndWeights(m::Vector{Float64}, V::AbstractMatrix; alpha=default_alpha, beta=default_beta, kappa=default_kappa)
    d = length(m)
    lambda = (d + kappa)*alpha^2 - d

    sigma_points = Vector{Vector{Float64}}(undef, 2*d+1)
    weights_m = Vector{Float64}(undef, 2*d+1)
    weights_c = Vector{Float64}(undef, 2*d+1)

    if isa(V, Diagonal)
        L = sqrt((d + lambda)*V) # Matrix square root
    else
        L = sqrt(Hermitian((d + lambda)*V))
    end

    sigma_points[1] = m
    weights_m[1] = lambda/(d + lambda)
    weights_c[1] = weights_m[1] + (1 - alpha^2 + beta)
    for i = 1:d
        sigma_points[2*i] = m + L[:,i]
        sigma_points[2*i+1] = m - L[:,i]
    end
    weights_m[2:end] .= 1/(2*(d + lambda))
    weights_c[2:end] .= 1/(2*(d + lambda))

    return (sigma_points, weights_m, weights_c)
end

"""
Return the statistics for the unscented approximation to the forward joint
"""
function unscentedStatistics(m::Float64, V::Float64, g::Function; alpha=default_alpha, beta=default_beta, kappa=default_kappa)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(m, V; alpha=alpha, beta=beta, kappa=kappa)

    g_sigma = g.(sigma_points)
    m_tilde = sum(weights_m.*g_sigma)
    V_tilde = sum(weights_c.*(g_sigma .- m_tilde).^2)
    C_tilde = sum(weights_c.*(sigma_points .- m).*(g_sigma .- m_tilde))

    return (m_tilde, V_tilde, C_tilde)
end

function unscentedStatistics(m::Vector{Float64}, V::AbstractMatrix, g::Function; alpha=default_alpha, beta=default_beta, kappa=default_kappa)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(m, V; alpha=alpha, beta=beta, kappa=kappa)
    d = length(m)

    g_sigma = g.(sigma_points)
    m_tilde = sum([weights_m[k+1]*g_sigma[k+1] for k=0:2*d])
    V_tilde = sum([weights_c[k+1]*(g_sigma[k+1] - m_tilde)*(g_sigma[k+1] - m_tilde)' for k=0:2*d])
    C_tilde = sum([weights_c[k+1]*(sigma_points[k+1] - m)*(g_sigma[k+1] - m_tilde)' for k=0:2*d])

    return (m_tilde, V_tilde, C_tilde)
end


#-------------
# Update Rules
#-------------

# Forward rule
function ruleSPNonlinearOutNG(msg_out::Nothing,
                              msg_in1::Message{F, V},
                              g::Function;
                              alpha::Float64=default_alpha) where {F<:Gaussian, V<:VariateType}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    (m_tilde, V_tilde, _) = unscentedStatistics(m_fw_in1, V_fw_in1, g; alpha=alpha)

    return Message(V, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

# Backward rule with given inverse
function ruleSPNonlinearIn1GG(msg_out::Message{F, V},
                              msg_in1::Nothing,
                              g::Function,
                              g_inv::Function;
                              alpha::Float64=default_alpha) where {F<:Gaussian, V<:VariateType}
    
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    (m_tilde, V_tilde, _) = unscentedStatistics(m_bw_out, V_bw_out, g_inv; alpha=alpha)

    return Message(V, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

# Backward rule with unknown inverse
function ruleSPNonlinearIn1GG(msg_out::Message{F1, V},
                              msg_in1::Message{F2, V},
                              g::Function;
                              alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian, V<:VariateType}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    W_fw_in1 = unsafePrecision(msg_in1.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)

    (m_tilde, V_tilde, C_tilde) = unscentedStatistics(m_fw_in1, V_fw_in1, g; alpha=alpha)

    # Update based on (Petersen et al. 2018; On Approximate Nonlinear Gaussian Message Passing on Factor Graphs)
    # Note, this implementation is not as efficient as Petersen et al. (2018), because we explicitly require the outbound messages
    C_tilde_inv = pinv(C_tilde)
    V_bw_in1 = V_fw_in1*C_tilde_inv'*(V_tilde + V_bw_out)*C_tilde_inv*V_fw_in1 - V_fw_in1
    m_bw_in1 = m_fw_in1 - (V_fw_in1 + V_bw_in1)*W_fw_in1*C_tilde*cholinv(V_tilde + V_bw_out)*(m_tilde - m_bw_out)

    return Message(V, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end


#---------------------------
# Custom inbounds collectors
#---------------------------

function collectSumProductNodeInbounds(node::Nonlinear, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    inbounds = Any[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if (node_interface == entry.interface == node.interfaces[2]) && (node.g_inv == nothing)
            # Collect the message inbound on the out edge if no inverse is available
            haskey(interface_to_schedule_entry, inbound_interface) || error("The nonlinear node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbounds, nothing)
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    # Push function (and inverse) to calling signature
    # These functions needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))
    if (entry.interface == node.interfaces[2]) && (node.g_inv != nothing)
        push!(inbounds, Dict{Symbol, Any}(:g_inv => node.g_inv,
                                          :keyword => false))
    end

    # Push spread parameter if manually defined
    if node.alpha != nothing
        push!(inbounds, Dict{Symbol, Any}(:alpha => node.alpha,
                                          :keyword => true))
    end

    return inbounds
end