export
ruleSPNonlinearSOutNM,
ruleSPNonlinearSIn1MN,
ruleSPNonlinearSOutNGX,
ruleSPNonlinearSInGX,
ruleMNonlinearSInGX,
prod!

const default_n_samples = 1000 # Default value for the number of samples


#----------------------
# Sampling Update Rules
#----------------------

function ruleSPNonlinearSOutNM(g::Function,
                               msg_out::Nothing,
                               msg_in1::Message{F, V};
                               n_samples=default_n_samples) where {F<:FactorFunction, V<:VariateType}

    samples = g.(sample(msg_in1.dist, n_samples))
    weights = ones(n_samples)/n_samples

    return Message(V, SampleList, s=samples, w=weights)
end

function ruleSPNonlinearSIn1MN(g::Function,
                               msg_out::Message{F, V},
                               msg_in1::Nothing;
                               n_samples=default_n_samples) where {F<:FactorFunction, V<:VariateType}

    return Message(V, Function, log_pdf = (z)->logPdf(msg_out.dist, g(z)))
end

function ruleSPNonlinearSOutNGX(g::Function,
                                msg_out::Nothing,
                                msgs_in::Vararg{Message{<:Gaussian, V}};
                                n_samples=default_n_samples) where V<:VariateType

    samples_in = [sample(msg_in.dist, n_samples) for msg_in in msgs_in]

    samples = g.(samples_in...)
    weights = ones(n_samples)/n_samples

    return Message(V, SampleList, s=samples, w=weights)
end

function ruleSPNonlinearSInGX(g::Function,
                              inx::Int64, # Index of inbound interface inx
                              msg_out::Message{<:Gaussian, V},
                              msgs_in::Vararg{Message{<:Gaussian, V}};
                              n_samples=default_n_samples) where V<:VariateType
    
    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    W_fw_in = cholinv(V_fw_in) # Convert to canonical statistics

    # Construct joint log-pdf function and gradient
    (log_joint, d_log_joint) = logJointPdfs(V, m_fw_in, W_fw_in, msg_out.dist, g, ds) # Overloaded on VariateType V

    # Compute joint belief on in's by gradient ascent
    m_in = gradientOptimization(log_joint, d_log_joint, m_fw_in, 0.01)
    V_in = cholinv(-ForwardDiff.jacobian(d_log_joint, m_in))

    # Marginalize joint belief on in's
    (m_inx, V_inx) = marginalizeGaussianMV(V, m_in, V_in, ds, inx) # Marginalization is overloaded on VariateType V
    W_inx = cholinv(V_inx) # Convert to canonical statistics
    xi_inx = W_inx*m_inx

    # Divide marginal on inx by forward message
    (xi_fw_inx, W_fw_inx) = unsafeWeightedMeanPrecision(msgs_in[inx].dist)
    xi_bw_inx = xi_inx - xi_fw_inx
    W_bw_inx = W_inx - W_fw_inx # Note: subtraction might lead to posdef inconsistencies

    return Message(V, GaussianWeightedMeanPrecision, xi=xi_bw_inx, w=W_bw_inx)
end

function ruleMNonlinearSInGX(g::Function,
                             msg_out::Message{<:Gaussian, V},
                             msgs_in::Vararg{Message{<:Gaussian, V}}) where V<:VariateType
    
    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    W_fw_in = cholinv(V_fw_in) # Convert to canonical statistics

    # Construct log-pdf function and gradient
    (log_joint, d_log_joint) = logJointPdfs(V, m_fw_in, W_fw_in, msg_out.dist, g, ds) # Overloaded on VariateType V

    # Compute joint marginal belief on in's by gradient ascent
    m_in = gradientOptimization(log_joint, d_log_joint, m_fw_in, 0.01)
    W_in = -ForwardDiff.jacobian(d_log_joint, m_in)

    return ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=m_in, w=W_in)
end


#---------------------------
# Custom product definitions
#---------------------------

function prod!(
    x::ProbabilityDistribution{V, Function},
    y::ProbabilityDistribution{V, Function},
    z::ProbabilityDistribution{V, Function}=ProbabilityDistribution(V, Function, log_pdf=(s)->s)) where V<:VariateType

    z.params[:log_pdf] = ( (s) -> x.params[:log_pdf](s) + y.params[:log_pdf](s) )

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate}, # Includes function distributions
    y::ProbabilityDistribution{Univariate, F},
    z::ProbabilityDistribution{Univariate, GaussianMeanPrecision}=ProbabilityDistribution(Univariate, GaussianMeanPrecision, m=0.0, w=1.0)) where {F<:Gaussian}

    # Optimize with gradient ascent
    log_joint(s) = logPdf(y,s) + logPdf(x,s)
    d_log_joint(s) = ForwardDiff.derivative(log_joint, s)
    m_initial = unsafeMean(y)
    
    mean = gradientOptimization(log_joint, d_log_joint, m_initial, 0.01)
    prec = -ForwardDiff.derivative(d_log_joint, mean)

    z.params[:m] = mean
    z.params[:w] = prec

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate}, # Includes function distributions
    y::ProbabilityDistribution{Multivariate, F},
    z::ProbabilityDistribution{Multivariate, GaussianMeanPrecision}=ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=[0.0], w=mat(1.0))) where {F<:Gaussian}

    # Optimize with gradient ascent
    log_joint(s) = logPdf(y,s) + logPdf(x,s)
    d_log_joint(s) = ForwardDiff.gradient(log_joint, s)
    m_initial = unsafeMean(y)

    mean = gradientOptimization(log_joint, d_log_joint, m_initial, 0.01)
    prec = -ForwardDiff.jacobian(d_log_joint, mean)

    z.params[:m] = mean
    z.params[:w] = prec

    return z
end


#--------------------------
# Custom inbound collectors
#--------------------------

function collectSumProductNodeInbounds(node::Nonlinear{Sampling}, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    inbounds = Any[]

    # Push function to calling signature; function needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))

    multi_in = (length(node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node
    inx = findfirst(isequal(entry.interface), node.interfaces) - 1 # Find number of inbound interface; 0 for outbound
    if multi_in && (inx > 0) # Multiple inbounds, and a backward message is required
        push!(inbounds, Dict{Symbol, Any}(:inx => inx, # Push inbound identifier to calling signature
                                          :keyword => false))
    end

    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if (node_interface == entry.interface != node.interfaces[1]) && multi_in
            # Collect the message inbound for backward rule
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

    if node.n_samples != nothing
        push!(inbounds, Dict{Symbol, Any}(:n_samples => node.n_samples,
                                          :keyword   => true))
    end

    return inbounds
end

function collectMarginalNodeInbounds(node::Nonlinear{Sampling}, entry::MarginalEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry
    inbound_cluster = entry.target # Entry target is a cluster

    inbounds = Any[]

    # Push function to calling signature; function needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))

    entry_pf = posteriorFactor(first(entry.target.edges))
    encountered_external_regions = Set{Region}()
    for node_interface in entry.target.node.interfaces
        current_region = region(inbound_cluster.node, node_interface.edge) # Note: edges that are not assigned to a posterior factor are assumed mean-field
        current_pf = posteriorFactor(node_interface.edge) # Returns an Edge if no posterior factor is assigned
        inbound_interface = ultimatePartner(node_interface)

        if (inbound_interface !== nothing) && isa(inbound_interface.node, Clamp)
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

    return inbounds
end


#---------------------------------
# Gradient optimization subroutine
#---------------------------------

function gradientOptimization(log_joint::Function, d_log_joint::Function, m_initial, step_size)
    dim_tot = length(m_initial)
    m_total = zeros(dim_tot)
    m_average = zeros(dim_tot)
    m_new = zeros(dim_tot)
    m_old = m_initial
    satisfied = false
    step_count = 0

    while !satisfied
        m_new = m_old .+ step_size.*d_log_joint(m_old)
        if log_joint(m_new) > log_joint(m_old)
            proposal_step_size = 10*step_size
            m_proposal = m_old .+ proposal_step_size.*d_log_joint(m_old)
            if log_joint(m_proposal) > log_joint(m_new)
                m_new = m_proposal
                step_size = proposal_step_size
            end
        else
            step_size = 0.1*step_size
            m_new = m_old .+ step_size.*d_log_joint(m_old)
        end
        step_count += 1
        m_total .+= m_old
        m_average = m_total ./ step_count
        if step_count > 10
            if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < dim_tot*0.1
                satisfied = true
            end
        end
        if step_count > dim_tot*250
            satisfied = true
        end
        m_old = m_new
    end

    return m_new
end


#--------
# Helpers
#--------

function logJointPdfs(V::Type{Multivariate}, m_fw_in::Vector, W_fw_in::AbstractMatrix, dist_out::ProbabilityDistribution, g::Function, ds::Vector{Int64})
    log_joint(x) = -0.5*sum(ds)*log(2pi) + 0.5*log(det(W_fw_in)) - 0.5*(x - m_fw_in)'*W_fw_in*(x - m_fw_in) + logPdf(dist_out, g(split(x, ds)...))
    d_log_joint(x) = ForwardDiff.gradient(log_joint, x)

    return (log_joint, d_log_joint)
end

function logJointPdfs(V::Type{Univariate}, m_fw_in::Vector, W_fw_in::AbstractMatrix, dist_out::ProbabilityDistribution, g::Function, ds::Vector{Int64})
    log_joint(x) = -0.5*sum(ds)*log(2pi) + 0.5*log(det(W_fw_in)) - 0.5*(x - m_fw_in)'*W_fw_in*(x - m_fw_in) + logPdf(dist_out, g(x...))
    d_log_joint(x) = ForwardDiff.gradient(log_joint, x)

    return (log_joint, d_log_joint)
end