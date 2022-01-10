export
ruleSPNonlinearCOutNM,
ruleSPNonlinearCIn1MN,
ruleSPNonlinearCOutNMX,
ruleSPNonlinearCInMX,
ruleSPNonlinearCInGX,
ruleMNonlinearCInMGX

const default_opt = ForgetDelayDescent(200.0, 0.6) # Default optimizer TODO: make variable
const default_n_iterations = 1000 # Default number of iterations for gradient descent # TODO: make variable


#-----------------------
# Conjugate Update Rules
#-----------------------

ruleSPNonlinearCOutNM(g::Function,
                      msg_out::Nothing,
                      msg_in1::Message;
                      dims::Any=nothing,
                      n_samples=default_n_samples) =
    ruleSPNonlinearSOutNM(g, nothing, msg_in1; dims=dims, n_samples=n_samples) # Reuse sampling update

function ruleSPNonlinearCIn1MN(g::Function,
                               msg_out::Message,
                               msg_in1::Message{F, V};
                               dims::Any=nothing,
                               n_samples=default_n_samples) where {F<:FactorNode, V<:VariateType}

    msg_s = ruleSPNonlinearSIn1MN(g, msg_out, nothing; dims=dims, n_samples=n_samples) # Returns Message{Function}
    η = naturalParams(msg_in1.dist)
    λ = renderCVI(msg_s.dist.params[:log_pdf], default_n_iterations, default_opt, η, msg_in1)

    return Message(standardDistribution(V, F, η=λ-η))
end

ruleSPNonlinearCOutNMX(g::Function,
                       msg_out::Nothing,
                       msgs_in::Vararg{Message};
                       dims::Any=nothing,
                       n_samples=default_n_samples) = 
    ruleSPNonlinearSOutNMX(g, nothing, msgs_in...; dims=dims, n_samples=n_samples)                                 

function ruleSPNonlinearCInGX(g::Function,
                              inx::Int64, # Index of inbound interface inx
                              msg_out::Message,
                              msgs_in::Vararg{Message{<:Gaussian}}; # Only Gaussian because of marginalization over inbounds
                              dims::Any=nothing,
                              n_samples=default_n_samples)

    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    msg_fw_in = Message(Multivariate, GaussianMeanVariance, m=m_fw_in, v=V_fw_in) # Joint forward message

    # log-pdf of joint backward message over inbounds
    log_pdf_s(z) = logPdf(msg_out.dist, g(split(z, ds)...))

    # Compute joint marginal belief
    η = naturalParams(msg_fw_in.dist)
    λ = renderCVI(log_pdf_s, default_n_iterations, default_opt, η, msg_fw_in)
    d_marg = standardDistribution(Multivariate, GaussianMeanVariance, η=λ)
    (m_in, V_in) = unsafeMeanCov(d_marg)
    
    # Marginalize joint belief on in's
    (m_inx, V_inx) = marginalizeGaussianMV(m_in, V_in, ds, inx)
    W_inx = cholinv(V_inx) # Convert to canonical statistics
    xi_inx = W_inx*m_inx

    # Divide marginal on inx by forward message
    (xi_fw_inx, W_fw_inx) = unsafeWeightedMeanPrecision(msgs_in[inx].dist)
    xi_bw_inx = xi_inx - xi_fw_inx
    W_bw_inx = W_inx - W_fw_inx # Note: subtraction might lead to posdef violations

    return Message(variateType(dims), GaussianWeightedMeanPrecision, xi=xi_bw_inx, w=W_bw_inx)
end

# Special case for two inputs with one PointMass (no inx required)
function ruleSPNonlinearCInMX(g::Function,
                              msg_out::Message,
                              msg_in1::Message{F, V},
                              msg_in2::Message{PointMass};
                              dims::Any=nothing,
                              n_samples=default_n_samples) where {F<:FactorNode, V<:VariateType}
    
    msg_s = ruleSPNonlinearSInMX(g, msg_out, nothing, msg_in2; dims=dims, n_samples=n_samples)
    η = naturalParams(msg_in1.dist)
    λ = renderCVI(msg_s.dist.params[:log_pdf], default_n_iterations, default_opt, η, msg_in1)

    return Message(standardDistribution(V, F, η=λ-η))
end

# Special case for two inputs with one PointMass (no inx required)
function ruleSPNonlinearCInMX(g::Function,
                              msg_out::Message,
                              msg_in1::Message{PointMass},
                              msg_in2::Message{F, V};
                              dims::Any=nothing,
                              n_samples=default_n_samples) where {F<:FactorNode, V<:VariateType}
    
    msg_s = ruleSPNonlinearSInMX(g, msg_out, msg_in1, nothing; dims=dims, n_samples=n_samples)
    η = naturalParams(msg_in2.dist)
    λ = renderCVI(msg_s.dist.params[:log_pdf], default_n_iterations, default_opt, η, msg_in2)

    return Message(standardDistribution(V, F, η=λ-η))
end

# Joint marginal belief over inbounds
function ruleMNonlinearCInMGX(g::Function,
                              msg_out::Message,
                              msgs_in::Vararg{Message{<:Gaussian}}) # Only Gaussian because of marginalization over inbounds
    
    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    msg_fw_in = Message(Multivariate, GaussianMeanVariance, m=m_fw_in, v=V_fw_in) # Joint forward message

    # log-pdf of joint backward message over inbounds
    log_pdf_s(z) = logPdf(msg_out.dist, g(split(z, ds)...))

    η = naturalParams(msg_fw_in.dist)
    λ = renderCVI(log_pdf_s, default_n_iterations, default_opt, η, msg_fw_in) # Natural statistics of marginal

    return standardDistribution(Multivariate, GaussianMeanVariance, η=λ)
end


#---------------------------
# Custom inbounds collectors
#---------------------------

# Conjugate approximation
function collectSumProductNodeInbounds(node::Nonlinear{Conjugate}, entry::ScheduleEntry)
    inbounds = Any[]

    # Push function to calling signature
    # This function needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))

    multi_in = isMultiIn(node) # Boolean to indicate a nonlinear node with multiple stochastic inbounds
    inx = findfirst(isequal(entry.interface), node.interfaces) - 1 # Find number of inbound interface; 0 for outbound
    
    if (inx > 0) && multi_in # Multi-inbound backward rule
        push!(inbounds, Dict{Symbol, Any}(:inx => inx, # Push inbound identifier
                                          :keyword => false))
    end

    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if (node_interface == entry.interface != node.interfaces[1]) && multi_in
            # Collect the breaker message for a backward rule with multiple inbounds
            haskey(interface_to_schedule_entry, inbound_interface) || error("The nonlinear node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbounds, nothing)
        elseif isClamped(inbound_interface)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    # Push custom arguments if defined
    if (node.dims !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:dims => node.dims[inx + 1],
                                          :keyword => true))
    end
    if (node.n_samples !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:n_samples => node.n_samples,
                                          :keyword   => true))
    end
    return inbounds
end


#-------------------------
# Optimization subroutines
#-------------------------

function renderCVI(logp_nc::Function,
                   n_its::Int,
                   opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent},
                   λ_init::Vector,
                   msg_in::Message{F, Univariate}) where F<:Gaussian

    η = deepcopy(naturalParams(msg_in.dist))
    λ = deepcopy(λ_init)

    df_m(z) = ForwardDiff.derivative(logp_nc, z)
    df_v(z) = 0.5*ForwardDiff.derivative(df_m, z)

    for i=1:n_its
        q = standardDistribution(Univariate, F, η=λ)
        z_s = sample(q)
        df_μ1 = df_m(z_s) - 2*df_v(z_s)*mean(q)
        df_μ2 = df_v(z_s)
        ∇f = [df_μ1, df_μ2]
        λ_old = deepcopy(λ)
        ∇ = λ .- η .- ∇f
        update!(opt, λ, ∇)
        if !isProper(standardDistribution(Univariate, F, η=λ))
            λ = λ_old
        end
    end

    return λ
end

function renderCVI(logp_nc::Function,
                   n_its::Int,
                   opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent},
                   λ_init::Vector,
                   msg_in::Message{F, Multivariate}) where F<:Gaussian

    η = deepcopy(naturalParams(msg_in.dist))
    λ = deepcopy(λ_init)

    df_m(z) = ForwardDiff.gradient(logp_nc, z)
    df_v(z) = 0.5*ForwardDiff.jacobian(df_m, z)

    for i=1:n_its
        q = standardDistribution(Multivariate, F, η=λ)
        z_s = sample(q)
        df_μ1 = df_m(z_s) - 2*df_v(z_s)*mean(q)
        df_μ2 = df_v(z_s)
        ∇f = [df_μ1; vec(df_μ2)]
        λ_old = deepcopy(λ)
        ∇ = λ .- η .- ∇f
        update!(opt, λ, ∇)
        if !isProper(standardDistribution(Multivariate, F, η=λ))
            λ = λ_old
        end
    end

    return λ
end

function renderCVI(logp_nc::Function,
                   n_its::Int,
                   opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent},
                   λ_init::Vector,
                   msg_in::Message{F, V}) where {F<:FactorNode, V<:VariateType}

    η = deepcopy(naturalParams(msg_in.dist)) # Concatenate natural parameters to vector
    λ = deepcopy(λ_init)

    A(η) = logNormalizer(msg_in.dist, η)
    gradA(η) = A'(η) # Zygote
    Fisher(η) = ForwardDiff.jacobian(gradA, η) # Zygote throws mutating array error
    for i=1:n_its
        q = standardDistribution(V, F, η=λ)
        z_s = sample(q)
        logq(λ) = logPdf(q, λ, z_s)
        ∇logq = logq'(λ)
        ∇f = Fisher(λ)\(logp_nc(z_s).*∇logq)
        λ_old = deepcopy(λ)
        ∇ = λ .- η .- ∇f
        update!(opt, λ, ∇)
        if !isProper(standardDistribution(V, F, η=λ))
            λ = λ_old
        end
    end

    return λ
end