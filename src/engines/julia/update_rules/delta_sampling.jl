export
ruleSPDeltaSOutNM,
ruleSPDeltaSIn1MN,
ruleSPDeltaSInGX,
ruleSPDeltaSOutNMX,
ruleSPDeltaSInMX,
ruleMDeltaSInMGX,
prod!

const default_n_samples = 1000 # Default value for the number of samples


#----------------------
# Sampling Update Rules
#----------------------

function ruleSPDeltaSOutNM(g::Function,
                           msg_out::Nothing,
                           msg_in1::Message; # Applies to any message except SampleList
                           dims::Any=nothing,
                           n_samples=default_n_samples)

    samples = g.(sample(msg_in1.dist, n_samples))
    weights = ones(n_samples)/n_samples

    return Message(variateType(dims), SampleList, s=samples, w=weights)
end

function ruleSPDeltaSOutNM(g::Function,
                           msg_out::Nothing,
                           msg_in1::Message{SampleList}; # Special case for SampleList
                           dims::Any=nothing,
                           n_samples=default_n_samples)

    samples = g.(msg_in1.dist.params[:s])
    weights = msg_in1.dist.params[:w]

    return Message(variateType(dims), SampleList, s=samples, w=weights)
end

function ruleSPDeltaSIn1MN(g::Function,
                           msg_out::Message,
                           msg_in1::Nothing;
                           dims::Any=nothing,
                           n_samples=default_n_samples)

    return Message(variateType(dims), Function, log_pdf = (z)->logPdf(msg_out.dist, g(z)))
end

function ruleSPDeltaSInGX(g::Function,
                          inx::Int64, # Index of inbound interface inx
                          msg_out::Message,
                          msgs_in::Vararg{Message{<:Gaussian}};
                          dims::Any=nothing,
                          n_samples=default_n_samples)

    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    W_fw_in = cholinv(V_fw_in) # Convert to canonical statistics

    # Construct joint log-pdf function and gradient
    (log_joint, d_log_joint) = logJointPdfs(m_fw_in, W_fw_in, msg_out.dist, g, ds)

    # Compute joint belief on in's by gradient ascent
    m_in = gradientOptimization(log_joint, d_log_joint, m_fw_in, 0.01)
    V_in = cholinv(-ForwardDiff.jacobian(d_log_joint, m_in))

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

function ruleSPDeltaSOutNMX(g::Function,
                            msg_out::Nothing,
                            msgs_in::Vararg{Message};
                            dims::Any=nothing,
                            n_samples=default_n_samples)
                                     
    samples_in = [sample(msg_in.dist, n_samples) for msg_in in msgs_in]
    samples = g.(samples_in...)
    weights = ones(n_samples)/n_samples

    return Message(variateType(dims), SampleList, s=samples, w=weights)
end

function ruleSPDeltaSInMX(g::Function,
                          inx::Int64, # Index of inbound interface inx
                          msg_out::Message,
                          msgs_in::Vararg{Message};
                          dims::Any=nothing,
                          n_samples=default_n_samples)

    arg_sample = (z) -> begin
        samples_in = []
        for i=1:length(msgs_in)
            if i==inx
                push!(samples_in, collect(Iterators.repeat([z], n_samples)))
            else
                push!(samples_in, sample(msgs_in[i].dist, n_samples))
            end
        end

        return samples_in
    end

    approximate_pdf(z) = sum(exp.(logPdf.([msg_out.dist],g.(arg_sample(z)...))))/n_samples

    return Message(variateType(dims), Function, log_pdf = (z)->log(approximate_pdf(z)))
end

# Special case for two inputs with one PointMass (no inx required)
function ruleSPDeltaSInMX(g::Function,
                          msg_out::Message,
                          msg_in1::Message{PointMass},
                          msg_in2::Nothing;
                          dims::Any=nothing,
                          n_samples=default_n_samples)

    m_in1 = msg_in1.dist.params[:m]

    return Message(variateType(dims), Function, log_pdf = (z)->logPdf(msg_out.dist, g(m_in1, z)))                          
end

# Special case for two inputs with one PointMass (no inx required)
function ruleSPDeltaSInMX(g::Function,
                          msg_out::Message,
                          msg_in1::Nothing,
                          msg_in2::Message{PointMass};
                          dims::Any=nothing,
                          n_samples=default_n_samples)

    m_in2 = msg_in2.dist.params[:m]

    return Message(variateType(dims), Function, log_pdf = (z)->logPdf(msg_out.dist, g(z, m_in2)))                          
end

function ruleMDeltaSInMGX(g::Function,
                          msg_out::Message,
                          msgs_in::Vararg{Message{<:Gaussian}})

    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    W_fw_in = cholinv(V_fw_in) # Convert to canonical statistics

    # Construct log-pdf function and gradient
    (log_joint, d_log_joint) = logJointPdfs(m_fw_in, W_fw_in, msg_out.dist, g, ds)

    # Compute joint marginal belief on in's by gradient ascent
    m_in = gradientOptimization(log_joint, d_log_joint, m_fw_in, 0.01)
    W_in = -ForwardDiff.jacobian(d_log_joint, m_in)

    return Distribution(Multivariate, Gaussian{Precision}, m=m_in, w=W_in)
end


#---------------------------
# Custom inbounds collectors
#---------------------------

# Sampling approximation
function collectSumProductNodeInbounds(node::Delta{Sampling}, entry::ScheduleEntry)
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
        if (node_interface == entry.interface != node.interfaces[1]) && multi_in
            # Collect the breaker message for a backward rule with multiple inbounds
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
    return inbounds
end


#---------------------------
# Custom product definitions
#---------------------------

function prod!(
    x::Distribution{V, Function},
    y::Distribution{V, Function}) where V<:VariateType # log-pdf for z cannot be predefined, because it cannot be overwritten

    return Distribution(V, Function, log_pdf=(s)->x.params[:log_pdf](s) + y.params[:log_pdf](s))
end

@symmetrical function prod!(
    x::Distribution{Univariate}, # Includes function distributions
    y::Distribution{Univariate, <:Gaussian},
    z::Distribution{Univariate, Gaussian{Precision}}=Distribution(Univariate, Gaussian{Precision}, m=0.0, w=1.0))

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
    x::Distribution{Multivariate}, # Includes function distributions
    y::Distribution{Multivariate, <:Gaussian},
    z::Distribution{Multivariate, Gaussian{Precision}}=Distribution(Multivariate, Gaussian{Precision}, m=[0.0], w=mat(1.0)))

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
    m_latests = if (dim_tot == 1) Queue{Float64}() else Queue{Vector}() end

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
        enqueue!(m_latests, m_old)
        if step_count > 10
            m_average = sum(x for x in m_latests)./10
            if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < dim_tot*0.1
                satisfied = true
            end
            dequeue!(m_latests);
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

function logJointPdfs(m_fw_in::Vector, W_fw_in::AbstractMatrix, dist_out::Distribution, g::Function, ds::Vector)
    log_joint(x) = -0.5*sum(intdim.(ds))*log(2pi) + 0.5*logdet(W_fw_in) - 0.5*(x - m_fw_in)'*W_fw_in*(x - m_fw_in) + logPdf(dist_out, g(split(x, ds)...))
    d_log_joint(x) = ForwardDiff.gradient(log_joint, x)

    return (log_joint, d_log_joint)
end