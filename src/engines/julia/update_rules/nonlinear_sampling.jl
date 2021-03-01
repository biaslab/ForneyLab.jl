export
ruleSPNonlinearSOutNM,
ruleSPNonlinearSIn1MN,
ruleSPNonlinearSOutNGX,
ruleSPNonlinearSInGX,
ruleSPNonlinearSOutNFactorX,
ruleSPNonlinearSInFactorX,
ruleMNonlinearSInGX,
prod!

const default_n_samples = 1000 # Default value for the number of samples


#----------------------
# Sampling Update Rules
#----------------------

function ruleSPNonlinearSOutNM(g::Function,
                               msg_out::Nothing,
                               msg_in1::Message{F, V};
                               n_samples=default_n_samples,
                               variate=V) where {F<:FactorFunction, V<:VariateType}

    samples = g.(sample(msg_in1.dist, n_samples))
    weights = ones(n_samples)/n_samples

    return Message(variate, SampleList, s=samples, w=weights)
end

function ruleSPNonlinearSIn1MN(g::Function,
                               msg_out::Message{F, V},
                               msg_in1::Nothing;
                               n_samples=default_n_samples,
                               variate=V) where {F<:FactorFunction, V<:VariateType}

    return Message(variate, Function, log_pdf = (z)->logPdf(msg_out.dist, g(z)))
end

function ruleSPNonlinearSOutNM(g::Function,
                               msg_out::Nothing,
                               msg_in1::Message{SampleList, V};
                               n_samples=default_n_samples,
                               variate=V) where {V<:VariateType}

    samples = g.(msg_in1.dist.params[:s])
    weights = msg_in1.dist.params[:w]

    return Message(variate, SampleList, s=samples, w=weights)
end

function msgSPNonlinearSOutNGX(g::Function,
    msg_out::Nothing,
    msgs_in::Vararg{Message{<:Gaussian, <:VariateType}};
    n_samples=default_n_samples,
    variate)

    samples_in = [sample(msg_in.dist, n_samples) for msg_in in msgs_in]

    samples = g.(samples_in...)
    weights = ones(n_samples)/n_samples

    return Message(variate, SampleList, s=samples, w=weights)
end

function ruleSPNonlinearSOutNGX(g::Function,
                                msg_out::Nothing,
                                msgs_in::Vararg{Message{<:Gaussian, <:VariateType}};
                                n_samples=default_n_samples,
                                variate)
    return msgSPNonlinearSOutNGX(g, msg_out, msgs_in..., n_samples=n_samples, variate=variate)
end

function ruleSPNonlinearSOutNGX(g::Function,
                                msg_out::Nothing,
                                msgs_in::Vararg{Message{<:Gaussian, V}};
                                n_samples=default_n_samples) where V<:VariateType
    return msgSPNonlinearSOutNGX(g, msg_out, msgs_in..., n_samples=n_samples, variate=V)
end

function msgSPNonlinearSInGX(g::Function,
                             inx::Int64, # Index of inbound interface inx
                             msg_out::Message{<:FactorFunction, <:VariateType},
                             msgs_in::Vararg{Message{<:Gaussian, <:VariateType}};
                             n_samples=default_n_samples,
                             variate)

    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    W_fw_in = cholinv(V_fw_in) # Convert to canonical statistics

    # Construct joint log-pdf function and gradient
    (log_joint, d_log_joint) = logJointPdfs(variate, m_fw_in, W_fw_in, msg_out.dist, g, ds) # Overloaded on VariateType V

    # Compute joint belief on in's by gradient ascent
    m_in = gradientOptimization(log_joint, d_log_joint, m_fw_in, 0.01)
    V_in = cholinv(-ForwardDiff.jacobian(d_log_joint, m_in))

    # Marginalize joint belief on in's
    (m_inx, V_inx) = marginalizeGaussianMV(variate, m_in, V_in, ds, inx) # Marginalization is overloaded on VariateType V
    W_inx = cholinv(V_inx) # Convert to canonical statistics
    xi_inx = W_inx*m_inx

    # Divide marginal on inx by forward message
    (xi_fw_inx, W_fw_inx) = unsafeWeightedMeanPrecision(msgs_in[inx].dist)
    xi_bw_inx = xi_inx - xi_fw_inx
    W_bw_inx = W_inx - W_fw_inx # Note: subtraction might lead to posdef inconsistencies

    return Message(variate, GaussianWeightedMeanPrecision, xi=xi_bw_inx, w=W_bw_inx)
end

function ruleSPNonlinearSInGX(g::Function,
                              inx::Int64, # Index of inbound interface inx
                              msg_out::Message{<:FactorFunction, V},
                              msgs_in::Vararg{Message{<:Gaussian, V}};
                              n_samples=default_n_samples) where V<:VariateType

    msgSPNonlinearSInGX(g, inx, msg_out, msgs_in..., n_samples=n_samples, variate=V)
end

function ruleSPNonlinearSInGX(g::Function,
                              inx::Int64, # Index of inbound interface inx
                              msg_out::Message{<:FactorFunction, <:VariateType},
                              msgs_in::Vararg{Message{<:Gaussian, <:VariateType}};
                              n_samples=default_n_samples,
                              variate)
    msgSPNonlinearSInGX(g, inx, msg_out, msgs_in..., n_samples=n_samples, variate=variate)
end

function msgSPNonlinearSOutNFactorX(g::Function,
                                    msg_out::Nothing,
                                    msgs_in::Vararg{Message{<:FactorNode}};
                                    n_samples=default_n_samples,
                                    variate)
    samples_in = [sample(msg_in.dist, n_samples) for msg_in in msgs_in]

    samples = g.(samples_in...)
    weights = ones(n_samples)/n_samples

    return Message(variate, SampleList, s=samples, w=weights)
end

function ruleSPNonlinearSOutNFactorX(g::Function,
                                     msg_out::Nothing,
                                     msgs_in::Vararg{Message{<:FactorNode, V}};
                                     n_samples=default_n_samples) where V<:VariateType
    msgSPNonlinearSOutNFactorX(g, msg_out, msgs_in..., n_samples=n_samples, variate=V)
end

function ruleSPNonlinearSOutNFactorX(g::Function,
                                     msg_out::Nothing,
                                     msgs_in::Vararg{Message{<:FactorNode}};
                                     n_samples=default_n_samples,
                                     variate)
    msgSPNonlinearSOutNFactorX(g, msg_out, msgs_in..., n_samples=n_samples, variate=variate)
end


function msgSPNonlinearSInFactorX(g::Function,
                                  inx::Int64, # Index of inbound interface inx
                                  msg_out::Message{<:FactorFunction},
                                  msgs_in::Vararg{Message{<:FactorNode}};
                                  n_samples=default_n_samples,
                                  variate)

    arg_sample = (z) -> begin
        samples_in = []
        for i=1:length(msgs_in)
            if i==inx
                push!(samples_in,collect(Iterators.repeat([ z ], n_samples)))
            else
                push!(samples_in,sample(msgs_in[i].dist, n_samples))
            end
        end

        return samples_in
    end

    approximate_pdf(z) = sum(exp.(logPdf.([msg_out.dist],g.(arg_sample(z)...))))/n_samples

    return Message(variate, Function, log_pdf = (z)->log(approximate_pdf(z)))
end

function ruleSPNonlinearSInFactorX(g::Function,
                                   inx::Int64, # Index of inbound interface inx
                                   msg_out::Message{<:FactorFunction, V},
                                   msgs_in::Vararg{Message{<:FactorNode, V}};
                                   n_samples=default_n_samples) where V<:VariateType

    msgSPNonlinearSInFactorX(g, inx, msg_out, msgs_in..., n_samples=n_samples, variate=V)

end

function ruleSPNonlinearSInFactorX(g::Function,
                                   inx::Int64, # Index of inbound interface inx
                                   msg_out::Message{<:FactorFunction},
                                   msgs_in::Vararg{Message{<:FactorNode}};
                                   n_samples=default_n_samples,
                                   variate)

    msgSPNonlinearSInFactorX(g, inx, msg_out, msgs_in..., n_samples=n_samples, variate=variate)
end

function ruleMNonlinearSInGX(g::Function,
                             msg_out::Message{<:FactorFunction, <:VariateType},
                             msgs_in::Vararg{Message{<:Gaussian, <:VariateType}})

    # Extract joint statistics of inbound messages
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Return arrays with individual means and covariances
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in) # Concatenate individual statistics into joint statistics
    W_fw_in = cholinv(V_fw_in) # Convert to canonical statistics

    # Construct log-pdf function and gradient
    (log_joint, d_log_joint) = logJointPdfs(Multivariate, m_fw_in, W_fw_in, msg_out.dist, g, ds) # Overloaded on VariateType V

    # Compute joint marginal belief on in's by gradient ascent
    m_in = gradientOptimization(log_joint, d_log_joint, m_fw_in, 0.01)
    W_in = -ForwardDiff.jacobian(d_log_joint, m_in)

    return ProbabilityDistribution(Multivariate, GaussianMeanPrecision, m=m_in, w=W_in)
end


#---------------------------
# Custom inbounds collectors
#---------------------------

# Unscented transform and extended approximation
function collectSumProductNodeInbounds(node::Nonlinear{Sampling}, entry::ScheduleEntry)
    inbounds = Any[]

    # Push function to calling signature
    # This function needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))

    multi_in = (length(node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node
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
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    # Push custom arguments if manually defined
    if (node.n_samples !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:n_samples => node.n_samples,
                                          :keyword   => true))
    end
    # Message on out interface
    if (inx == 0) && (node.out_variate !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:variate => node.out_variate,
                                          :keyword   => true))
    end
    # Message on in interface
    if (inx > 0) && (node.in_variates !== nothing)
        push!(inbounds, Dict{Symbol, Any}(:variate => node.in_variates[inx],
                                          :keyword   => true))
    end
    return inbounds
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
    m_latests = if dim_tot == 1 Queue{Float64}() else Queue{Vector}() end

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

function logJointPdfs(V::Type{Multivariate}, m_fw_in::Vector, W_fw_in::AbstractMatrix, dist_out::ProbabilityDistribution, g::Function, ds::Vector{Int64})
    log_joint(x) = -0.5*sum(ds)*log(2pi) + 0.5*logdet(W_fw_in) - 0.5*(x - m_fw_in)'*W_fw_in*(x - m_fw_in) + logPdf(dist_out, g(split(x, ds)...))
    d_log_joint(x) = ForwardDiff.gradient(log_joint, x)

    return (log_joint, d_log_joint)
end

function logJointPdfs(V::Type{Univariate}, m_fw_in::Vector, W_fw_in::AbstractMatrix, dist_out::ProbabilityDistribution, g::Function, ds::Vector{Int64})
    log_joint(x) = -0.5*sum(ds)*log(2pi) + 0.5*logdet(W_fw_in) - 0.5*(x - m_fw_in)'*W_fw_in*(x - m_fw_in) + logPdf(dist_out, g(x...))
    d_log_joint(x) = ForwardDiff.gradient(log_joint, x)

    return (log_joint, d_log_joint)
end
