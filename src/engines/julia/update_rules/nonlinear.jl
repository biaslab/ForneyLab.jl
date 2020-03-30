export
ruleSPNonlinearUTOutNG,
ruleSPNonlinearUTIn1GG,
ruleSPNonlinearISInMN,
ruleSPNonlinearISOutNG,
ruleSPNonlinearLInMN,
ruleSPNonlinearLOutNG,
prod!


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

function ruleSPNonlinearUTOutNG(msg_out::Nothing,
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

function ruleSPNonlinearUTOutNG(msg_out::Nothing,
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

function ruleSPNonlinearUTIn1GG(msg_out::Message{F, Univariate},
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

function ruleSPNonlinearUTIn1GG(msg_out::Message{F, Multivariate},
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

function ruleSPNonlinearUTIn1GG(msg_out::Message{F1, Univariate},
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

function ruleSPNonlinearUTIn1GG(msg_out::Message{F1, Multivariate},
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

function collectSumProductNodeInbounds(node::Nonlinear{Unscented}, entry::ScheduleEntry)
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

function ruleSPNonlinearISInMN(msg_out::Message{F, Univariate}, msg_in1::Nothing, g::Function) where {F<:SoftFactor}
    Message(Univariate, Function, log_pdf=(z) -> logPdf(msg_out.dist, g(z)), ApproximationType="NonlinearIS")
end

function ruleSPNonlinearISOutNG(msg_out::Nothing, msg_in1::Message{F, Univariate}, g::Function) where {F<:Gaussian}
    # The forward message is parameterized by a SampleList
    dist_in1 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)

    sample_list = g.(dist_in1.params[:m] .+ sqrt(dist_in1.params[:v]).*randn(1000))

    weight_list = ones(1000)/1000

    Message(Univariate, SampleList, s=sample_list, w=weight_list)
end

function ruleSPNonlinearLInMN(msg_out::Message{F, Univariate}, msg_in1::Nothing, g::Function) where {F<:SoftFactor}
    try
        ForwardDiff.derivative(g, 0)
        return Message(Univariate, Function, log_pdf=(z) -> logPdf(msg_out.dist, g(z)), ApproximationType="NonlinearL")
    catch
        return Message(Multivariate, Function, log_pdf=(z) -> logPdf(msg_out.dist, g(z)), ApproximationType="NonlinearL")
    end
end

function ruleSPNonlinearLOutNG(msg_out::Nothing, msg_in1::Message{F, Univariate}, g::Function) where {F<:Gaussian}
    # The forward message is parameterized by a SampleList
    dist_in1 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)

    sample_list = g.(dist_in1.params[:m] .+ sqrt(dist_in1.params[:v]).*randn(1000))

    weight_list = ones(1000)/1000

    return Message(Univariate, SampleList, s=sample_list, w=weight_list)

    # if length(g(dist_in1.params[:m])) == 1
    #     sample_list = g.(dist_in1.params[:m] .+ sqrt(dist_in1.params[:v]).*randn(1000))
    #
    #     return Message(Univariate, SampleList, s=sample_list)
    # else
    #     sample_list = g.(dist_in1.params[:m] .+ sqrt(dist_in1.params[:v]).*randn(1000))
    #
    #     return Message(Multivariate, SampleList, s=sample_list)
    # end
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, F},
    z::ProbabilityDistribution{Univariate, GaussianMeanVariance}=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) where {F<:Gaussian}

    if x.params[:ApproximationType] == "NonlinearIS"
        # The product of a log-pdf and Gaussian distribution is computed by importance sampling
        y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
        samples = y.params[:m] .+ sqrt(y.params[:v]).*randn(1000)

        p = exp.((x.params[:log_pdf]).(samples))
        Z = sum(p)
        mean = sum(p./Z.*samples)
        var = sum(p./Z.*(samples .- mean).^2)

        z.params[:m] = mean
        z.params[:v] = var
    elseif x.params[:ApproximationType] == "NonlinearL"
        # The product of a log-pdf and Gaussian distribution is approximated by Laplace method
        y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
        log_joint(s) = logPdf(y,s) + x.params[:log_pdf](s)
        #Optimization with gradient ascent
        d_log_joint(s) = ForwardDiff.derivative(log_joint, s)
        m_old = y.params[:m] #initial point
        step_size = 0.01 #initial step size
        satisfied = 0
        step_count = 0
        m_total = 0.0
        m_average = 0.0
        m_new = 0.0
        while satisfied == 0
            m_new = m_old + step_size*d_log_joint(m_old)
            if log_joint(m_new) > log_joint(m_old)
                proposal_step_size = 10*step_size
                m_proposal = m_old + proposal_step_size*d_log_joint(m_old)
                if log_joint(m_proposal) > log_joint(m_new)
                    m_new = m_proposal
                    step_size = proposal_step_size
                end
            else
                step_size = 0.1*step_size
                m_new = m_old + step_size*d_log_joint(m_old)
            end
            step_count += 1
            m_total += m_old
            m_average = m_total / step_count
            if step_count > 10
                if abs((m_new-m_average)/m_average) < 0.1
                    satisfied = 1
                end
            elseif step_count > 250
                satisfied = 1
            end
            m_old = m_new
        end
        mean = m_new
        var = - 1.0 / ForwardDiff.derivative(d_log_joint, mean)

        z.params[:m] = mean
        z.params[:v] = var
    elseif x.params[:ApproximationType] == "BivariateL"
        z.params[:m] = x.params[:m]
        z.params[:v] = x.params[:v]
    end

    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate, Function},
    y::ProbabilityDistribution{Multivariate, F},
    z::ProbabilityDistribution{Multivariate, GaussianMeanVariance}=ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=[0.0], v=mat(1.0))) where {F<:Gaussian}

    if x.params[:ApproximationType] == "NonlinearL"
        # The product of a log-pdf and Gaussian distribution is approximated by Laplace method
        y = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, y)
        dim = dims(y)
        log_joint(s) = logPdf(y,s) + x.params[:log_pdf](s)
        #Optimization with gradient ascent
        d_log_joint(s) = ForwardDiff.gradient(log_joint, s)
        m_old = y.params[:m] #initial point
        step_size = 0.01 #initial step size
        satisfied = 0
        step_count = 0
        m_total = zeros(dim)
        m_average = zeros(dim)
        m_new = zeros(dim)
        while satisfied == 0
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
                if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < dim*0.1
                    satisfied = 1
                end
            elseif step_count > dim*250
                satisfied = 1
            end
            m_old = m_new
        end
        mean = m_new
        var = inv(- 1.0 .* ForwardDiff.jacobian(d_log_joint, mean))

        z.params[:m] = mean
        z.params[:v] = var
    elseif x.params[:ApproximationType] == "BivariateL"
        z.params[:m] = x.params[:m]
        z.params[:v] = x.params[:v]
    end

    return z
end

# Think more carefully about the prod of function messages
function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, Function},
    z::ProbabilityDistribution{Univariate, Function}=ProbabilityDistribution(Univariate, Function, log_pdf=(s)->s, ApproximationType="NonlinearIS"))

    z.params[:log_pdf] = ((s) -> x.params[:log_pdf](s) + y.params[:log_pdf](s))
    if x.params[:ApproximationType] == y.params[:ApproximationType]
        z.params[:ApproximationType] = x.params[:ApproximationType]
    else
        z.params[:ApproximationType] = "NonlinearIS"
    end

    return z
end

function prod!(
    x::ProbabilityDistribution{Multivariate, Function},
    y::ProbabilityDistribution{Multivariate, Function},
    z::ProbabilityDistribution{Multivariate, Function}=ProbabilityDistribution(Multivariate, Function, log_pdf=(s)->s, ApproximationType="NonlinearL"))

    z.params[:log_pdf] = ((s) -> x.params[:log_pdf](s) + y.params[:log_pdf](s))
    if x.params[:ApproximationType] == y.params[:ApproximationType]
        z.params[:ApproximationType] = x.params[:ApproximationType]
    else
        z.params[:ApproximationType] = "NonlinearL"
    end

    return z
end

#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::Nonlinear{ImportanceSampling}, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    inbounds = Any[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
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

    return inbounds
end

function collectSumProductNodeInbounds(node::Nonlinear{Laplace}, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    inbounds = Any[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
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

    return inbounds
end
