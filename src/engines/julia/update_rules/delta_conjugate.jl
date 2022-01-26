# Inplementation based on Khan et al. (2017), "Conjugate-computation variational inference: 
# Converting variational inference in non-conjugate models to inferences in conjugate models",
# and Akbayrak et al. (2021), "Extended Variational Message Passing for Automated Approximate
# Bayesian Inference"

export
ruleSPDeltaCOutNM,
ruleSPDeltaCIn1MN,
ruleSPDeltaCOutNMX,
ruleSPDeltaCInMX,
ruleSPDeltaCInGX,
ruleMDeltaCInMGX

const default_optimizer = ForgetDelayDescent(200.0, 0.6) # Default optimizer
const default_n_iterations = 1000 # Default number of iterations for gradient descent


#-----------------------
# Conjugate Update Rules
#-----------------------

ruleSPDeltaCOutNM(g::Function,
                  msg_out::Nothing,
                  msg_in1::Message;
                  dims::Any=nothing,
                  n_samples=default_n_samples,
                  n_iterations=default_n_iterations,
                  optimizer=default_optimizer) =
    ruleSPDeltaSOutNM(g, nothing, msg_in1; dims=dims, n_samples=n_samples) # Reuse sampling update

function ruleSPDeltaCIn1MN(g::Function,
                           msg_out::Message,
                           msg_in1::Message{F, V};
                           dims::Any=nothing,
                           n_samples=default_n_samples,
                           n_iterations=default_n_iterations,
                           optimizer=default_optimizer) where {F<:FactorNode, V<:VariateType}

    msg_s = ruleSPDeltaSIn1MN(g, msg_out, nothing; dims=dims, n_samples=n_samples) # Returns Message{Function}
    η = naturalParams(msg_in1.dist)
    λ = renderCVI(msg_s.dist.params[:log_pdf], n_iterations, optimizer, η, msg_in1)

    return Message(standardDistribution(V, F, η=λ-η))
end

ruleSPDeltaCOutNMX(g::Function,
                   msg_out::Nothing,
                   msgs_in::Vararg{Message};
                   dims::Any=nothing,
                   n_samples=default_n_samples,
                   n_iterations=default_n_iterations,
                   optimizer=default_optimizer) =
     ruleSPDeltaSOutNMX(g, nothing, msgs_in...; dims=dims, n_samples=n_samples)                                 

function ruleSPDeltaCInGX(g::Function,
                          inx::Int64, # Index of inbound interface inx
                          msg_out::Message,
                          msgs_in::Vararg{Message{<:Gaussian}}; # Only Gaussian because of marginalization over inbounds
                          dims::Any=nothing,
                          n_samples=default_n_samples,
                          n_iterations=default_n_iterations,
                          optimizer=default_optimizer)                                      

    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    msg_fw_in = Message(Multivariate, Gaussian{Moments}, m=m_fw_in, v=V_fw_in) # Joint forward message

    # log-pdf of joint backward message over inbounds
    log_pdf_s(z) = logPdf(msg_out.dist, g(split(z, ds)...))

    # Compute joint marginal belief
    η = naturalParams(msg_fw_in.dist)
    λ = renderCVI(log_pdf_s, n_iterations, optimizer, η, msg_fw_in)
    d_marg = standardDistribution(Multivariate, Gaussian{Moments}, η=λ)
    (m_in, V_in) = unsafeMeanCov(d_marg)
    
    # Marginalize joint belief on in's
    (m_inx, V_inx) = marginalizeGaussianMV(m_in, V_in, ds, inx)
    W_inx = cholinv(V_inx) # Convert to canonical statistics
    xi_inx = W_inx*m_inx

    # Divide marginal on inx by forward message
    (xi_fw_inx, W_fw_inx) = unsafeWeightedMeanPrecision(msgs_in[inx].dist)
    xi_bw_inx = xi_inx - xi_fw_inx
    W_bw_inx = W_inx - W_fw_inx # Note: subtraction might lead to posdef violations

    return Message(variateType(dims), Gaussian{Canonical}, xi=xi_bw_inx, w=W_bw_inx)
end

# Special case for two inputs with one PointMass (no inx required)
function ruleSPDeltaCInMX(g::Function,
                          msg_out::Message,
                          msg_in1::Message{F, V},
                          msg_in2::Message{PointMass};
                          dims::Any=nothing,
                          n_samples=default_n_samples,
                          n_iterations=default_n_iterations,
                          optimizer=default_optimizer) where {F<:FactorNode, V<:VariateType}
    
    msg_s = ruleSPDeltaSInMX(g, msg_out, nothing, msg_in2; dims=dims, n_samples=n_samples)
    η = naturalParams(msg_in1.dist)
    λ = renderCVI(msg_s.dist.params[:log_pdf], n_iterations, optimizer, η, msg_in1)

    return Message(standardDistribution(V, F, η=λ-η))
end

# Special case for two inputs with one PointMass (no inx required)
function ruleSPDeltaCInMX(g::Function,
                          msg_out::Message,
                          msg_in1::Message{PointMass},
                          msg_in2::Message{F, V};
                          dims::Any=nothing,
                          n_samples=default_n_samples,
                          n_iterations=default_n_iterations,
                          optimizer=default_optimizer) where {F<:FactorNode, V<:VariateType}
    
    msg_s = ruleSPDeltaSInMX(g, msg_out, msg_in1, nothing; dims=dims, n_samples=n_samples)
    η = naturalParams(msg_in2.dist)
    λ = renderCVI(msg_s.dist.params[:log_pdf], n_iterations, optimizer, η, msg_in2)

    return Message(standardDistribution(V, F, η=λ-η))
end

# Joint marginal belief over inbounds
function ruleMDeltaCInMGX(g::Function,
                          msg_out::Message,
                          msgs_in::Vararg{Message{<:Gaussian}}) # Only Gaussian because of marginalization over inbounds
    
    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    msg_fw_in = Message(Multivariate, Gaussian{Moments}, m=m_fw_in, v=V_fw_in) # Joint forward message

    # log-pdf of joint backward message over inbounds
    log_pdf_s(z) = logPdf(msg_out.dist, g(split(z, ds)...))

    η = naturalParams(msg_fw_in.dist)
    λ = renderCVI(log_pdf_s, default_n_iterations, default_optimizer, η, msg_fw_in) # Natural statistics of marginal

    return standardDistribution(Multivariate, Gaussian{Moments}, η=λ)
end


#---------------------------
# Custom inbounds collectors
#---------------------------

# Conjugate approximation
function collectSumProductNodeInbounds(node::Delta{Conjugate}, entry::ScheduleEntry)
    inbounds = Any[]

    # Push function to calling signature
    # This function needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))

    multi_in = isMultiIn(node) # Boolean to indicate a Delta node with multiple stochastic inbounds
    inx = findfirst(isequal(entry.interface), node.interfaces) - 1 # Find number of inbound interface; 0 for outbound
    
    if (inx > 0) && multi_in # Multi-inbound backward rule
        push!(inbounds, Dict{Symbol, Any}(:inx => inx, # Push inbound identifier
                                          :keyword => false))
    end

    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if (node_interface == entry.interface != node.interfaces[1])
            # Collect the breaker message for a backward rule
            haskey(interface_to_schedule_entry, inbound_interface) || error("The Delta node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
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
    if (node.n_iterations !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:n_iterations => node.n_iterations,
                                          :keyword      => true))
    end
    if (node.optimizer !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:optimizer => node.optimizer,
                                          :keyword   => true))
    end
    return inbounds
end


#-------------------------
# Optimization subroutines
#-------------------------

function renderCVI(log_μ_bw::Function,
                   n_iterations::Int,
                   optimizer::Any,
                   λ_0::Vector,
                   μ_fw::Message{F, V}) where {F<:FactorNode, V<:VariateType}

    # Intialize natural parameters of forward message
    η = naturalParams(μ_fw.dist)                   

    # Initialize Fisher information matrix
    A = λ -> logNormalizer(V, F, η=λ)
    Fisher = λ -> ForwardDiff.hessian(A, λ)

    # Initialize q marginal
    λ_i = deepcopy(λ_0)
    q_i = standardDistribution(V, F, η=λ_i)

    for i=1:n_iterations
        # Store previous results for possible reset
        q_i_min = deepcopy(q_i)
        λ_i_min = deepcopy(λ_i)

        # Given the current sample, define natural gradient of q
        s_q_i = sample(q_i)
        log_q = λ -> logPdf(V, F, s_q_i, η=λ)
        ∇log_q = λ -> ForwardDiff.gradient(log_q, λ)

        # Compute current free energy gradient and update natural statistics
        ∇log_μ_bw_i = log_μ_bw(s_q_i)*cholinv(Fisher(λ_i))*∇log_q(λ_i) # Natural gradient of backward message
        ∇F_i = λ_i - η - ∇log_μ_bw_i # Natural gradient of free energy
        λ_i -= apply!(optimizer, λ_i, ∇F_i) # Update λ_i

        # Update q_i
        q_i = standardDistribution(V, F, η=λ_i)
        if !isProper(q_i) # Result is improper; reset statistics
            q_i = q_i_min
            λ_i = λ_i_min
        end
    end

    return λ_i
end

# Gaussian result that avoids Fisher information matrix construction
function renderCVI(log_μ_bw::Function,
                   n_iterations::Int,
                   optimizer::Any,
                   λ_0::Vector,
                   μ_fw::Message{F, V}) where {F<:Gaussian, V<:VariateType}

    # Intialize natural parameters of forward message
    η = naturalParams(μ_fw.dist)                   

    # Intialize gradients/derivatives of Gaussian moments
    if V == Univariate
        ∇m = s -> ForwardDiff.derivative(log_μ_bw, s)
        ∇v = s -> 0.5*ForwardDiff.derivative(∇m, s)
    else
        ∇m = s -> ForwardDiff.gradient(log_μ_bw, s)
        ∇v = s -> 0.5*ForwardDiff.jacobian(∇m, s)
    end

    # Initialize q marginal
    λ_i = deepcopy(λ_0)
    q_i = standardDistribution(V, F, η=λ_i)

    for i=1:n_iterations
        # Store previous results for possible reset
        q_i_min = deepcopy(q_i)
        λ_i_min = deepcopy(λ_i)

        # Given the current sample, define natural gradient of q
        m_q_i = unsafeMean(q_i)
        s_q_i = sample(q_i)
        ∇λ_i_1 = ∇m(s_q_i) - 2*∇v(s_q_i)*m_q_i
        ∇λ_i_2 = ∇v(s_q_i)
        
        # Compute current free energy gradient and update natural statistics
        ∇log_μ_bw_i = vcat(∇λ_i_1, vec(∇λ_i_2))
        ∇F_i = λ_i - η - ∇log_μ_bw_i # Natural gradient of free energy
        λ_i -= apply!(optimizer, λ_i, ∇F_i) # Update λ_i

        # Update q_i
        q_i = standardDistribution(V, F, η=λ_i)
        if !isProper(q_i) # Result is improper; reset statistics
            q_i = q_i_min
            λ_i = λ_i_min
        end
    end

    return λ_i
end