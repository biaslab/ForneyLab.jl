export
ruleSPRGMPInFN,
ruleSPRGMPOutNG,
prod!

function ruleSPRGMPInFN(msg_out::Message,msg_z::Nothing,g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64)
    log_grad_g(z) = log(abs(ForwardDiff.derivative(g,z)))
    Message(Univariate, Function, log_pdf=(z)-> log_grad_g(z) + logPdf(msg_out.dist,g(z)), eta_a=lr, eta_b=lr, num_epoch=num_epochs1, num_epoch2=num_epochs2)
end

function ruleSPRGMPOutNG(msg_out::Nothing,msg_z::Message{F1, V},g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64) where {F1<:Gaussian, V<:Univariate}
    #z = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_z)
    z_m = msg_z.dist.params[:m]
    z_v = msg_z.dist.params[:v]
    Message(V, RGMP_dist, m=z_m, v=z_v, f=g)
end

#New prod rules for RGMP
#differentiable reparameterization function for Gaussian
function g_func(epsilon, a, b)
    return a + b*epsilon
end

function gaussian_avi1(log_fA, log_fB, mu_prior, v_prior, eta_a, eta_b, num_epoch, num_epoch2)
    # We use the prior parameters as the initial point of the variational parameters
    # Initial parameters can be different
    a_t = mu_prior
    b_t = sqrt(v_prior)
    # There are two epochs, the first one is where we update the variational parameters
    # The second epoch is for monte carlo gradients which can be set to 1 for fully stochastic optimization
    log_f(x) = log_fA(x) + log_fB(x)
    for epoch=1:num_epoch
        sum_a = 0
        sum_b = 0
        for epoch2=1:num_epoch2
            epsilon = randn()
            logf_a(a) = log_f(g_func(epsilon, a, b_t))
            logf_b(b) = log_f(g_func(epsilon, a_t, b))
            grad_logf_a(a) = ForwardDiff.derivative(logf_a,a)
            grad_logf_b(b) = ForwardDiff.derivative(logf_b,b)
            a_grad = grad_logf_a(a_t)
            #a_second_grad = ForwardDiff.derivative(grad_logf_a,a_t)
            #a_grad = -(1/a_second_grad)*a_grad
            b_grad = grad_logf_b(b_t) + 1.0/b_t
            #b_second_grad = ForwardDiff.derivative(grad_logf_b,b_t)
            #b_grad = -(1/b_second_grad)*b_grad
            sum_a += a_grad
            sum_b += b_grad
        end
        a_t = a_t + eta_a*sum_a/num_epoch2
        b_t = b_t + eta_b*sum_b/num_epoch2
    end
    return a_t, b_t^2
end

function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, GaussianMeanVariance}=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) where {F2<:Gaussian}

    log_FB(s) = x.params[:log_pdf](s)
    #pdf of prior which refers to the upcoming message from the previous time dependent variable
    #fA(mu,v,s) = Distributions.pdf(Distributions.Normal(mu, v),s)
    y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
    y_m = y.params[:m]
    y_v = y.params[:v]
    #Insert the parameters of previous message into prior function
    #log_FA(s) = log(fA(y_m,y_v,s))
    log_FA(s) = logPdf(y,s)
    #approximate the posterior with Gaussian distribution family
    eta_a = x.params[:eta_a]
    eta_b = x.params[:eta_b]
    num_epoch = x.params[:num_epoch]
    num_epoch2 = x.params[:num_epoch2]
    mu_zt, v_zt = gaussian_avi1(log_FA, log_FB, y_m, y_v, eta_a, eta_b, num_epoch, num_epoch2)

    z.params[:m] = mu_zt
    z.params[:v] = v_zt

    return z
end

function prod!(
    y::ProbabilityDistribution{Univariate, F1},
    x::ProbabilityDistribution{Univariate, Function},
    z::ProbabilityDistribution{Univariate, GaussianMeanVariance}=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) where {F1<:Gaussian}

    #log pdf of likelihood function Po(exp(.))
    log_FB(s) = x.params[:log_pdf](s)
    #pdf of prior which refers to the upcoming message from the previous time dependent variable
    #fA(mu,v,s) = Distributions.pdf(Distributions.Normal(mu, v),s)
    y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
    y_m = y.params[:m]
    y_v = y.params[:v]
    #Insert the parameters of previous message into prior function
    #log_FA(s) = log(fA(y_m,y_v,s))
    log_FA(s) = logPdf(y,s)
    eta_a = x.params[:eta_a]
    eta_b = x.params[:eta_b]
    num_epoch = x.params[:num_epoch]
    num_epoch2 = x.params[:num_epoch2]
    mu_zt, v_zt = gaussian_avi1(log_FA, log_FB, y_m, y_v, eta_a, eta_b, num_epoch, num_epoch2)

    z.params[:m] = mu_zt
    z.params[:v] = v_zt

    return z
end

function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, Function},
    z::ProbabilityDistribution{Univariate, Function}=ProbabilityDistribution(Univariate, Function, log_pdf=(z)-> z, eta_a=0.01, eta_b=0.01, num_epoch=1, num_epoch2=1000))

    log_FA = x.params[:log_pdf]
    log_FB = y.params[:log_pdf]
    #approximate the posterior with Gaussian distribution family
    eta_a = x.params[:eta_a]
    eta_b = x.params[:eta_b]
    num_epoch = x.params[:num_epoch]
    num_epoch2 = x.params[:num_epoch2]

    z.params[:log_pdf] = (s -> log_FA(s) + log_FB(s))
    z.params[:eta_a] = eta_a
    z.params[:eta_b] = eta_b
    z.params[:num_epoch] = num_epoch
    z.params[:num_epoch2] = num_epoch2

    return z
end

#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::RGMP, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
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

    push!(inbound_messages, "$(node.g)")
    push!(inbound_messages, "$(node.num_epochs1)")
    push!(inbound_messages, "$(node.num_epochs2)")
    push!(inbound_messages, "$(node.lr)")

    return inbound_messages
end
