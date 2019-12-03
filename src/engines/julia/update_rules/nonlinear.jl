export
ruleSPNonlinearOutNG,
ruleSPNonlinearIn1GG

# Determine the default value for the spread parameter
const default_alpha = 1e-3

"""
Return the sigma points and weights for a Gaussian distribution
"""
function sigmaPointsAndWeights(dist::ProbabilityDistribution{Univariate, F}, alpha::Float64) where F<:Gaussian
    (m_x, V_x) = unsafeMeanCov(dist)

    kappa = 0
    beta = 2
    lambda = (1 + kappa)*alpha^2 - 1

    sigma_points = Vector{Float64}(undef, 3)
    weights_m = Vector{Float64}(undef, 3)
    weights_c = Vector{Float64}(undef, 3)

    l = sqrt((1 + lambda)*V_x)

    sigma_points[1] = m_x
    sigma_points[2] = m_x + l
    sigma_points[3] = m_x - l
    weights_m[1] = lambda/(1 + lambda)
    weights_m[2] = weights_m[3] = 1/(2*(1 + lambda))
    weights_c[1] = weights_m[1] + (1 - alpha^2 + beta)
    weights_c[2] = weights_c[3] = 1/(2*(1 + lambda))

    return (sigma_points, weights_m, weights_c)
end

function sigmaPointsAndWeights(dist::ProbabilityDistribution{Multivariate, F}, alpha::Float64) where F<:Gaussian
    d = dims(dist)
    (m_x, V_x) = unsafeMeanCov(dist)

    kappa = 0
    beta = 2
    lambda = (d + kappa)*alpha^2 - d

    sigma_points = Vector{Vector{Float64}}(undef, 2*d+1)
    weights_m = Vector{Float64}(undef, 2*d+1)
    weights_c = Vector{Float64}(undef, 2*d+1)

    if isa(V_x, Diagonal)
        L = sqrt((d + lambda)*V_x) # Matrix square root
    else
        L = sqrt(Hermitian((d + lambda)*V_x))
    end

    sigma_points[1] = m_x
    weights_m[1] = lambda/(d + lambda)
    weights_c[1] = weights_m[1] + (1 - alpha^2 + beta)
    for i = 1:d
        sigma_points[2*i] = m_x + L[:,i]
        sigma_points[2*i+1] = m_x - L[:,i]
    end
    weights_m[2:end] .= 1/(2*(d + lambda))
    weights_c[2:end] .= 1/(2*(d + lambda))

    return (sigma_points, weights_m, weights_c)
end

function ruleSPNonlinearOutNG(msg_out::Nothing,
                              msg_in1::Message{F, Univariate},
                              g::Function;
                              alpha::Float64=default_alpha) where F<:Gaussian

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist, alpha)

    # Unscented approximation
    g_sigma = g.(sigma_points)
    m_fw_out = sum(weights_m.*g_sigma)
    V_fw_out = sum(weights_c.*(g_sigma .- m_fw_out).^2)

    return Message(Univariate, GaussianMeanVariance, m=m_fw_out, v=V_fw_out)
end

function ruleSPNonlinearOutNG(msg_out::Nothing,
                              msg_in1::Message{F, Multivariate},
                              g::Function;
                              alpha::Float64=default_alpha) where F<:Gaussian
    d = dims(msg_in1.dist)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist, alpha)

    # Unscented approximation
    g_sigma = g.(sigma_points)
    m_fw_out = sum([weights_m[k+1]*g_sigma[k+1] for k=0:2*d])
    V_fw_out = sum([weights_c[k+1]*(g_sigma[k+1] - m_fw_out)*(g_sigma[k+1] - m_fw_out)' for k=0:2*d])

    return Message(Multivariate, GaussianMeanVariance, m=m_fw_out, v=V_fw_out)
end

function ruleSPNonlinearIn1GG(msg_out::Message{F, Univariate},
                              msg_in1::Nothing,
                              g::Function,
                              g_inv::Function;
                              alpha::Float64=default_alpha) where F<:Gaussian

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_out.dist, alpha)

    # Unscented approximation
    g_inv_sigma = g_inv.(sigma_points)
    m_bw_in1 = sum(weights_m.*g_inv_sigma)
    V_bw_in1 = sum(weights_c.*(g_inv_sigma .- m_bw_in1).^2)

    return Message(Univariate, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end

function ruleSPNonlinearIn1GG(msg_out::Message{F, Multivariate},
                              msg_in1::Nothing,
                              g::Function,
                              g_inv::Function;
                              alpha::Float64=default_alpha) where F<:Gaussian
    d = dims(msg_out.dist)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_out.dist, alpha)

    # Unscented approximation
    g_inv_sigma = g_inv.(sigma_points)
    m_bw_in1 = sum([weights_m[k+1]*g_inv_sigma[k+1] for k=0:2*d])
    V_bw_in1 = sum([weights_c[k+1]*(g_inv_sigma[k+1] - m_bw_in1)*(g_inv_sigma[k+1] - m_bw_in1)' for k=0:2*d])

    return Message(Multivariate, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end

function ruleSPNonlinearIn1GG(msg_out::Message{F1, Univariate},
                              msg_in1::Message{F2, Univariate},
                              g::Function;
                              alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist, alpha)

    # Unscented approximations
    g_sigma = g.(sigma_points)
    m_fw_out = sum(weights_m.*g_sigma)
    V_fw_out = sum(weights_c.*(g_sigma .- m_fw_out).^2)
    C_fw = sum(weights_c.*(sigma_points .- m_fw_in1).*(g_sigma .- m_fw_out))

    # Update based on (Petersen et al. 2018; On Approximate Nonlinear Gaussian Message Passing on Factor Graphs)
    C_fw_inv = 1/C_fw
    V_bw_in1 = V_fw_in1^2*C_fw_inv^2*(V_fw_out + V_bw_out) - V_fw_in1
    m_bw_in1 = m_fw_in1 - (V_fw_in1 + V_bw_in1)*C_fw*(m_fw_out - m_bw_out)/(V_fw_in1*(V_fw_out + V_bw_out))

    return Message(Univariate, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end

function ruleSPNonlinearIn1GG(msg_out::Message{F1, Multivariate},
                              msg_in1::Message{F2, Multivariate},
                              g::Function;
                              alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian}
    d_in1 = dims(msg_in1.dist)

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    W_fw_in1 = unsafePrecision(msg_in1.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist, alpha)

    # Unscented approximations
    g_sigma = g.(sigma_points)
    m_fw_out = sum([weights_m[k+1]*g_sigma[k+1] for k=0:2*d_in1])
    V_fw_out = sum([weights_c[k+1]*(g_sigma[k+1] - m_fw_out)*(g_sigma[k+1] - m_fw_out)' for k=0:2*d_in1])
    C_fw = sum([weights_c[k+1]*(sigma_points[k+1] - m_fw_in1)*(g_sigma[k+1] - m_fw_out)' for k=0:2*d_in1])

    # Update based on (Petersen et al. 2018; On Approximate Nonlinear Gaussian Message Passing on Factor Graphs)
    # Note, this implementation is not as efficient as Petersen et al. (2018), because we explicitly require the outbound messages
    C_fw_inv = pinv(C_fw)
    V_bw_in1 = V_fw_in1*C_fw_inv'*(V_fw_out + V_bw_out)*C_fw_inv*V_fw_in1 - V_fw_in1
    m_bw_in1 = m_fw_in1 - (V_fw_in1 + V_bw_in1)*W_fw_in1*C_fw*cholinv(V_fw_out + V_bw_out)*(m_fw_out - m_bw_out)

    return Message(Multivariate, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::Nonlinear, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if (node_interface == entry.interface == node.interfaces[2]) && (node.g_inv == nothing)
            # Collect the message inbound on the out edge if no inverse is available
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

    # Push function (and inverse) to calling signature
    # These functions needs to be defined in the scope of the user
    push!(inbound_messages, "$(node.g)")
    if (entry.interface == node.interfaces[2]) && (node.g_inv != nothing)
        push!(inbound_messages, "$(node.g_inv)")
    end

    # Push spread parameter if manually defined
    if node.alpha != nothing
        push!(inbound_messages, "alpha=$(node.alpha)")
    end

    return inbound_messages
end
