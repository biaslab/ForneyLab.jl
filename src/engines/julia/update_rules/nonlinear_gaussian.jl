export
ruleSPNonlinearGaussianOutNGP,
ruleSPNonlinearGaussianIn1GGP,
ruleMNonlinearGaussianGGD

const default_alpha = 1e-3 # Default value for the spread parameter
const default_beta = 2.0
const default_kappa = 0.0

"""
Return the sigma points and weights for a Gaussian distribution
"""
function sigmaPointsAndWeights(dist::ProbabilityDistribution{Univariate, F}; alpha=default_alpha, beta=default_beta, kappa=default_kappa) where F<:Gaussian
    (m, V) = unsafeMeanCov(dist)
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

function sigmaPointsAndWeights(dist::ProbabilityDistribution{Multivariate, F}; alpha=default_alpha, beta=default_beta, kappa=default_kappa) where F<:Gaussian
    d = dims(dist)
    (m, V) = unsafeMeanCov(dist)
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

function ruleSPNonlinearGaussianOutNGP(msg_out::Nothing,
                                       msg_in1::Message{F, Univariate},
                                       msg_var::Message{PointMass, Univariate},
                                       g::Function;
                                       alpha::Float64=default_alpha) where F<:Gaussian

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist; alpha=alpha)

    # Unscented approximation
    g_sigma = g.(sigma_points)
    m_tilde = sum(weights_m.*g_sigma)
    V_tilde = sum(weights_c.*(g_sigma .- m_tilde).^2) + msg_var.dist.params[:m] # Include additive noise

    return Message(Univariate, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

function ruleSPNonlinearGaussianOutNGP(msg_out::Nothing,
                                       msg_in1::Message{F, Multivariate},
                                       msg_var::Message{PointMass, MatrixVariate},
                                       g::Function;
                                       alpha::Float64=default_alpha) where F<:Gaussian
    d = dims(msg_in1.dist)
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist; alpha=alpha)

    # Unscented approximation
    g_sigma = g.(sigma_points)
    m_tilde = sum([weights_m[k+1]*g_sigma[k+1] for k=0:2*d])
    V_tilde = sum([weights_c[k+1]*(g_sigma[k+1] - m_tilde)*(g_sigma[k+1] - m_tilde)' for k=0:2*d]) + msg_var.dist.params[:m] # Include additive noise

    return Message(Multivariate, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

# Backward rule with given inverse
function ruleSPNonlinearGaussianIn1GGP(msg_out::Message{F, Univariate},
                                       msg_in1::Nothing,
                                       msg_var::Message{PointMass, Univariate},
                                       g::Function,
                                       g_inv::Function;
                                       alpha::Float64=default_alpha) where F<:Gaussian
    
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    noisy_dist_out = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=m_bw_out, v=V_bw_out+msg_var.dist.params[:m])

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(noisy_dist_out; alpha=alpha)

    # Unscented approximation
    g_inv_sigma = g_inv.(sigma_points)
    m_tilde = sum(weights_m.*g_inv_sigma)
    V_tilde = sum(weights_c.*(g_inv_sigma .- m_tilde).^2)

    return Message(Univariate, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

# Backward rule with given inverse
function ruleSPNonlinearGaussianIn1GGP(msg_out::Message{F, Multivariate},
                                       msg_in1::Nothing,
                                       msg_var::Message{PointMass, MatrixVariate},
                                       g::Function,
                                       g_inv::Function;
                                       alpha::Float64=default_alpha) where F<:Gaussian
    d = dims(msg_out.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    noisy_dist_out = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_bw_out, v=V_bw_out+msg_var.dist.params[:m])
    
    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(noisy_dist_out; alpha=alpha)

    # Unscented approximation
    g_inv_sigma = g_inv.(sigma_points)
    m_tilde = sum([weights_m[k+1]*g_inv_sigma[k+1] for k=0:2*d])
    V_tilde = sum([weights_c[k+1]*(g_inv_sigma[k+1] - m_tilde)*(g_inv_sigma[k+1] - m_tilde)' for k=0:2*d])

    return Message(Multivariate, GaussianMeanVariance, m=m_tilde, v=V_tilde)
end

function ruleSPNonlinearGaussianIn1GGP(msg_out::Message{F1, Univariate},
                                       msg_in1::Message{F2, Univariate},
                                       msg_var::Message{PointMass, Univariate},
                                       g::Function;
                                       alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian}

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist; alpha=alpha)

    # Unscented approximations
    g_sigma = g.(sigma_points)
    m_tilde = sum(weights_m.*g_sigma)
    V_tilde = sum(weights_c.*(g_sigma .- m_tilde).^2) + msg_var.dist.params[:m]
    C_tilde = sum(weights_c.*(sigma_points .- m_fw_in1).*(g_sigma .- m_tilde))

    # Update based on (Petersen et al. 2018; On Approximate NonlinearGaussian Gaussian Message Passing on Factor Graphs)
    C_tilde_inv = 1/C_tilde
    V_bw_in1 = V_fw_in1^2*C_tilde_inv^2*(V_tilde + V_bw_out) - V_fw_in1
    m_bw_in1 = m_fw_in1 - (V_fw_in1 + V_bw_in1)*C_tilde*(m_tilde - m_bw_out)/(V_fw_in1*(V_tilde + V_bw_out))

    return Message(Univariate, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end

# Backward rule with unknown inverse
function ruleSPNonlinearGaussianIn1GGP(msg_out::Message{F1, Multivariate},
                                       msg_in1::Message{F2, Multivariate},
                                       msg_var::Message{PointMass, MatrixVariate},
                                       g::Function;
                                       alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian}
    d_in1 = dims(msg_in1.dist)

    (m_fw_in1, V_fw_in1) = unsafeMeanCov(msg_in1.dist)
    W_fw_in1 = unsafePrecision(msg_in1.dist)
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist; alpha=alpha)

    # Unscented approximations
    g_sigma = g.(sigma_points)
    m_tilde = sum([weights_m[k+1]*g_sigma[k+1] for k=0:2*d_in1])
    V_tilde = sum([weights_c[k+1]*(g_sigma[k+1] - m_tilde)*(g_sigma[k+1] - m_tilde)' for k=0:2*d_in1]) + msg_var.dist.params[:m]
    C_tilde = sum([weights_c[k+1]*(sigma_points[k+1] - m_fw_in1)*(g_sigma[k+1] - m_tilde)' for k=0:2*d_in1])

    # Update based on (Petersen et al. 2018; On Approximate NonlinearGaussian Gaussian Message Passing on Factor Graphs)
    # Note, this implementation is not as efficient as Petersen et al. (2018), because we explicitly require the outbound messages
    C_tilde_inv = pinv(C_tilde)
    V_bw_in1 = V_fw_in1*C_tilde_inv'*(V_tilde + V_bw_out)*C_tilde_inv*V_fw_in1 - V_fw_in1
    m_bw_in1 = m_fw_in1 - (V_fw_in1 + V_bw_in1)*W_fw_in1*C_tilde*cholinv(V_tilde + V_bw_out)*(m_tilde - m_bw_out)

    return Message(Multivariate, GaussianMeanVariance, m=m_bw_in1, v=V_bw_in1)
end

# Backward rule with unknown inverse
function ruleMNonlinearGaussianGGD(msg_out::Message{F1, Univariate},
                                   msg_in1::Message{F2, Univariate},
                                   dist_var::ProbabilityDistribution{Univariate},
                                   g::Function;
                                   alpha::Float64=default_alpha) where {F1<:Gaussian, F2<:Gaussian}

    (xi_fw_in1, W_fw_in1) = unsafeWeightedMeanPrecision(msg_in1.dist)
    m_fw_in1 = unsafeMean(msg_in1.dist)
    (xi_bw_out, W_bw_out) = unsafeWeightedMeanPrecision(msg_out.dist)

    (sigma_points, weights_m, weights_c) = sigmaPointsAndWeights(msg_in1.dist; alpha=alpha)

    # Unscented approximations
    g_sigma = g.(sigma_points)
    m_tilde = sum(weights_m.*g_sigma)
    V_tilde = sum(weights_c.*(g_sigma .- m_tilde).^2) + unsafeMean(dist_var)
    C_tilde = sum(weights_c.*(sigma_points .- m_fw_in1).*(g_sigma .- m_tilde))
    
    # Compute joint statistics
    z = m_tilde - C_tilde'*xi_fw_in1
    G_tilde = 1/(V_tilde - C_tilde'*W_fw_in1*C_tilde)
    K = W_fw_in1*C_tilde*G_tilde
    xi_q_out = G_tilde*z + xi_bw_out
    xi_q_in1 = xi_fw_in1 - K*z
    W_q_out = G_tilde + W_bw_out
    W_q_in1 = W_fw_in1 + K*C_tilde'*W_fw_in1

    return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[xi_q_out; xi_q_in1], w=[W_q_out -K'; -K W_q_in1])
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::NonlinearGaussian, entry::ScheduleEntry)
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

function collectMarginalNodeInbounds(node::NonlinearGaussian, entry::MarginalEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry
    inbound_cluster = entry.target # Entry target is a cluster

    inbounds = Any[]
    entry_pf = posteriorFactor(first(entry.target.edges))
    encountered_external_regions = Set{Region}()
    for node_interface in entry.target.node.interfaces
        current_region = region(inbound_cluster.node, node_interface.edge) # Note: edges that are not assigned to a posterior factor are assumed mean-field 
        current_pf = posteriorFactor(node_interface.edge) # Returns an Edge if no posterior factor is assigned
        inbound_interface = ultimatePartner(node_interface)

        if (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Edge is clamped, hard-code marginal of constant node
            push!(inbounds, assembleClamp!(copy(inbound_interface.node), ProbabilityDistribution)) # Copy Clamp before assembly to prevent overwriting dist_or_msg field
        elseif (current_pf === entry_pf)
            # Edge is internal, collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif !(current_region in encountered_external_regions)
            # Edge is external and region is not yet encountered, collect marginal from marginal dictionary
            push!(inbounds, target_to_marginal_entry[current_region])
            push!(encountered_external_regions, current_region) # Register current region with encountered external regions
        end
    end

    # Push function to calling signature
    # These functions needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))

    # Push spread parameter if manually defined
    if node.alpha != nothing
        push!(inbounds, Dict{Symbol, Any}(:alpha => node.alpha,
                                          :keyword => true))
    end

    return inbounds
end
