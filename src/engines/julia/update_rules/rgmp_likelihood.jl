export
ruleSPRGMP_likelihoodIn1PN


#differentiable reparameterization function for Gaussian
function g_func(epsilon, a, b)
    return a + b*epsilon
end

#approximating exact marginal posterior with gaussian q and auto-vi
function gaussian_avi(log_fA, log_fB, mu_prior, v_prior, eta_a, eta_b, num_epoch, num_epoch2)
    # We use the prior parameters as the initial point of the variational parameters
    # Initial parameters can be different
    a_t = mu_prior
    b_t = sqrt(v_prior)
    # There are two epochs, the first one is where we update the variational parameters
    # The second epoch is for monte carlo gradients which can be set to 1 for fully stochastic optimization
    for epoch=1:num_epoch
        sum_a = 0
        sum_b = 0
        for epoch2=1:num_epoch2
            epsilon = randn()
            post_sample = g_func(epsilon, a_t, b_t)
            logfA_derivative = ForwardDiff.derivative(log_fA, post_sample)
            logfB_derivative = ForwardDiff.derivative(log_fB, post_sample)
            logf_derivative = logfA_derivative + logfB_derivative
            g_a(a) = g_func(epsilon, a, b_t)
            a_grad = logf_derivative * ForwardDiff.derivative(g_a, a_t)
            g_b(b) = g_func(epsilon, a_t, b)
            b_grad = logf_derivative * ForwardDiff.derivative(g_b, b_t) + 1.0/b_t
            sum_a += a_grad
            sum_b += b_grad
        end
        a_t = a_t + eta_a*sum_a/num_epoch2
        b_t = b_t + eta_b*sum_b/num_epoch2
    end
    return a_t, b_t^2
end

function gaussian_avi_mv_uv(log_f, mu_prior, v_prior, eta_a, eta_b, num_epoch, num_epoch2)
    # We use the prior parameters as the initial point of the variational parameters
    # Initial parameters can be different
    a_t = mu_prior
    d = length(a_t) #dimensionalty for the incoming MV Gaussian
    b_t = cholesky(v_prior).L #Lower triangular cholesky decomposition component
    # There are two epochs, the first one is where we update the variational parameters
    # The second epoch is for monte carlo gradients which can be set to 1 for fully stochastic optimization
    for epoch=1:num_epoch
        sum_a = zeros(d)
        sum_b = zeros(d,d)
        for epoch2=1:num_epoch2
            epsilon = randn(d)
            logf_a(a) = log_f(g_func(epsilon, a, b_t))
            logf_b(b) = log_f(g_func(epsilon, a_t, b))
            a_grad = ForwardDiff.gradient(logf_a,a_t)
            b_grad = ForwardDiff.gradient(logf_b,b_t) + Diagonal(1.0 ./ b_t)
            sum_a += a_grad
            sum_b += b_grad
        end
        a_t = a_t + eta_a*sum_a/num_epoch2
        b_t = b_t + eta_b*sum_b/num_epoch2
    end
    return a_t, b_t*transpose(b_t)
end

#mean and std. parameterization
univareateGaussianLogPdf(mu,s,x) = -log(sqrt(2*pi)) -log(s) - ((x-mu)^2)/(2*s^2)
#mean and covariance parameterization
multivariateGaussianLogPdf(mu,V,d,x) = -0.5*d*log(2*pi) -0.5*log(det(V)) -0.5*transpose(x-mu)*inv(V)*(x-mu)


function ruleSPRGMP_likelihoodIn1PN(msg_out::Message{PointMass, Univariate}, msg_in1::Message{GaussianMeanVariance, Univariate}, g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64)

    obs = msg_out.dist.params[:m]
    # We have an observation so insert observation value into log-likelhood function
    log_FB(x) = g(obs,x)
    z = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)
    z_m = z.params[:m]
    z_v = z.params[:v]
    #pdf of prior which refers to the upcoming message from the previous time dependent variable
    #Insert the parameters of previous message into prior function
    log_FA(x) = univareateGaussianLogPdf(z_m,sqrt(z_v),x)
    #approximate the posterior with Gaussian distribution family
    mu_zt, v_zt = gaussian_avi(log_FA, log_FB, z_m, z_v, lr, lr, num_epochs1, num_epochs2)
    v_message = 1/(1/v_zt - 1/z_v)
    m_message = v_message*(mu_zt/v_zt - z_m/z_v)
    Message(Univariate, GaussianMeanVariance, m=m_message, v=v_message)
end

function ruleSPRGMP_likelihoodIn1PN(msg_out::Message{PointMass, Multivariate}, msg_in1::Message{GaussianMeanVariance, Multivariate}, g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64)

    obs = msg_out.dist.params[:m]
    # We have an observation so insert observation value into log-likelhood function
    log_FB(x) = g(obs,x)
    z = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in1.dist)
    z_m = z.params[:m] #Supposed to be d element array
    z_v = z.params[:v]
    d = length(z_m) #dimensionalty for the incoming MV Gaussian
    #pdf of prior which refers to the upcoming message from the previous time dependent variable
    #Insert the parameters of previous message into prior function
    log_FA(x) = multivariateGaussianLogPdf(z_m,z_v,d,x)
    # Joint distribution function
    log_F(x) = log_FA(x) + log_FB(x)
    #approximate the posterior with Gaussian distribution family
    try
        #theoretically prior should be incorporated in the function to be optimized
        #but due to the design of ForneyLab, the message is approximated by dividing
        #the approximate posterior to the Gaussian prior which may cause losing the
        #positive definity constraint of covariance matrix
        mu_zt, v_zt = gaussian_avi_mv_uv(log_F, z_m, z_v, lr, lr, num_epochs1, num_epochs2)
        inv_v_zt, inv_z_v = inv(v_zt), inv(z_v)
        v_message = Hermitian(inv(inv_v_zt-inv_z_v))
        cholesky(v_message) #check if the covariance matrix is positive
        m_message = v_message*(inv_v_zt*mu_zt - inv_z_v*z_m)
        Message(Multivariate, GaussianMeanVariance, m=m_message, v=v_message)
    catch
        #if the division leads to non-positive definite covariance matrix then go
        #with approximation through likelihood function!
        mu_zt, v_zt = gaussian_avi_mv_uv(log_FB, z_m, z_v, lr, lr, num_epochs1, num_epochs2)
        Message(Multivariate, GaussianMeanVariance, m=mu_zt, v=v_zt)
    end
end

function ruleSPRGMP_likelihoodIn1PN(msg_out::Message{PointMass, Univariate}, msg_in1::Message{GaussianMeanVariance, Multivariate}, g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64)

    obs = msg_out.dist.params[:m]
    # We have an observation so insert observation value into log-likelhood function
    log_FB(x) = g(obs,x)
    z = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in1.dist)
    z_m = z.params[:m] #Supposed to be d element array
    z_v = z.params[:v]
    d = length(z_m) #dimensionalty for the incoming MV Gaussian
    #pdf of prior which refers to the upcoming message from the previous time dependent variable
    #Insert the parameters of previous message into prior function
    log_FA(x) = multivariateGaussianLogPdf(z_m,z_v,d,x)
    # Joint distribution function
    log_F(x) = log_FA(x) + log_FB(x)
    #approximate the posterior with Gaussian distribution family
    try
        #theoretically prior should be incorporated in the function to be optimized
        #but due to the design of ForneyLab, the message is approximated by dividing
        #the approximate posterior to the Gaussian prior which may cause losing the
        #positive definity constraint of covariance matrix
        mu_zt, v_zt = gaussian_avi_mv_uv(log_F, z_m, z_v, lr, lr, num_epochs1, num_epochs2)
        inv_v_zt, inv_z_v = inv(v_zt), inv(z_v)
        v_message = Hermitian(inv(inv_v_zt-inv_z_v))
        cholesky(v_message) #check if the covariance matrix is positive
        m_message = v_message*(inv_v_zt*mu_zt - inv_z_v*z_m)
        Message(Multivariate, GaussianMeanVariance, m=m_message, v=v_message)
    catch
        #if the division leads to non-positive definite covariance matrix then go
        #with approximation through likelihood function!
        mu_zt, v_zt = gaussian_avi_mv_uv(log_FB, z_m, z_v, lr, lr, num_epochs1, num_epochs2)
        Message(Multivariate, GaussianMeanVariance, m=mu_zt, v=v_zt)
    end
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::RGMP_likelihood, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == node.interfaces[2]
            # Always collect the message inbound on the in1 edge,
            # because it is used to obtain the approximation point
            haskey(interface_to_msg_idx, inbound_interface) || error("The nonlinear node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        elseif node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbound_messages, "nothing")
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbound_messages, ForneyLab.messageString(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        end
    end

    push!(inbound_messages, "$(node.g)")
    push!(inbound_messages, "$(node.num_epochs1)")
    push!(inbound_messages, "$(node.num_epochs2)")
    push!(inbound_messages, "$(node.lr)")

    return inbound_messages
end
