export
ruleSPBivariateLIn1MNG,
ruleSPBivariateLIn2MGN,
ruleSPBivariateLOutNGG,
ruleMBivariateLOutNGG

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F1, Univariate}, msg_in2::Message{F2, Univariate}, g::Function, status::Dict) where {F1<:Gaussian,F2<:Gaussian}
    # The forward message is parameterized by a SampleList
    dist_in1 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)
    dist_in2 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in2.dist)

    samples1 = dist_in1.params[:m] .+ sqrt(dist_in1.params[:v]).*randn(1000)
    samples2 = dist_in2.params[:m] .+ sqrt(dist_in2.params[:v]).*randn(1000)

    sample_list = g.(samples1, samples2)
    weight_list = ones(1000)/1000

    if length(sample_list[1]) == 1
        return Message(Univariate, SampleList, s=sample_list, w=weight_list)
    else
        return Message(Multivariate, SampleList, s=sample_list, w=weight_list)
    end
end

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F1, Multivariate}, msg_in2::Message{F2, Univariate}, g::Function, status::Dict) where {F1<:Gaussian,F2<:Gaussian}
    # The forward message is parameterized by a SampleList
    dist_in1 = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in1.dist)
    dist_in2 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in2.dist)

    C1L = cholesky(dist_in1.params[:v]).L
    dim = dims(dist_in1)
    samples1 = []
    for j=1:1000
        sample = dist_in1.params[:m] + C1L*randn(dim)
        push!(samples1,sample)
    end
    samples2 = dist_in2.params[:m] .+ sqrt(dist_in2.params[:v]).*randn(1000)

    sample_list = g.(samples1, samples2)
    weight_list = ones(1000)/1000

    if length(sample_list[1]) == 1
        return Message(Univariate, SampleList, s=sample_list, w=weight_list)
    else
        return Message(Multivariate, SampleList, s=sample_list, w=weight_list)
    end
end

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F1, Univariate}, msg_in2::Message{F2, Multivariate}, g::Function, status::Dict) where {F1<:Gaussian,F2<:Gaussian}
    # The forward message is parameterized by a SampleList
    dist_in1 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)
    dist_in2 = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in2.dist)

    samples1 = dist_in1.params[:m] .+ sqrt(dist_in1.params[:v]).*randn(1000)

    C2L = cholesky(dist_in2.params[:v]).L
    dim = dims(dist_in2)
    samples2 = []
    for j=1:1000
        sample = dist_in2.params[:m] + C2L*randn(dim)
        push!(samples2,sample)
    end

    sample_list = g.(samples1, samples2)
    weight_list = ones(1000)/1000

    if length(sample_list[1]) == 1
        return Message(Univariate, SampleList, s=sample_list, w=weight_list)
    else
        return Message(Multivariate, SampleList, s=sample_list, w=weight_list)
    end
end

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F1, Multivariate}, msg_in2::Message{F2, Multivariate}, g::Function, status::Dict) where {F1<:Gaussian,F2<:Gaussian}
    # The forward message is parameterized by a SampleList
    dist_in1 = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in1.dist)
    dist_in2 = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in2.dist)

    C1L = cholesky(dist_in1.params[:v]).L
    dim1 = dims(dist_in1)
    samples1 = []

    C2L = cholesky(dist_in2.params[:v]).L
    dim2 = dims(dist_in2)
    samples2 = []

    for j=1:1000
        sample1 = dist_in1.params[:m] + C1L*randn(dim1)
        sample2 = dist_in2.params[:m] + C2L*randn(dim2)
        push!(samples1,sample1)
        push!(samples2,sample2)
    end

    sample_list = g.(samples1, samples2)
    weight_list = ones(1000)/1000

    if length(sample_list[1]) == 1
        return Message(Univariate, SampleList, s=sample_list, w=weight_list)
    else
        return Message(Multivariate, SampleList, s=sample_list, w=weight_list)
    end
end

function approxMessageBivariate(m_prior::Number,v_prior::Number,m_post::Number,v_post::Number)

    if abs(v_prior-v_post) < 1e-5
        v_message = 1e-5
    else
        v_message = v_prior*v_post/(v_prior-v_post)
    end
    m_message = (m_post*(v_prior+v_message) - m_prior*v_message)/v_prior
    return Message(Univariate, GaussianMeanVariance, m=m_message, v=v_message)
end

function approxMessageBivariate(m_prior::Array,v_prior,m_post::Array,v_post)

    w_prior, w_post = inv(v_prior+2e-5*diageye(length(m_prior))), inv(v_post+1e-5*diageye(length(m_prior)))
    w_message = w_post - w_prior
    xi_message = (w_prior+w_message)*m_post - w_prior*m_prior
    return Message(Multivariate, GaussianWeightedMeanPrecision, xi=xi_message, w=w_message)
end


function ruleSPBivariateLIn1MNG(msg_out::Message{Fout, Vout}, msg_in1::Message{F1, V1}, msg_in2::Message{F2, V2}, g::Function, status::Dict) where {Fout<:SoftFactor, Vout<:VariateType, F1<:Gaussian, V1<:VariateType, F2<:Gaussian, V2<:VariateType}

    if status[:updated]
        status[:updated] = false
        return status[:message]
    else
        dist_in1 = convert(ProbabilityDistribution{V1, GaussianMeanVariance}, msg_in1.dist)
        dist_in2 = convert(ProbabilityDistribution{V2, GaussianMeanVariance}, msg_in2.dist)

        m_concat = [dist_in1.params[:m];dist_in2.params[:m]]

        dim1 = dims(dist_in1)
        dim2 = dims(dist_in2)
        dim_tot = dim1 + dim2
        v_concat = zeros(dim_tot,dim_tot)

        if dim1 == 1
            v_concat[dim1,dim1] = dist_in1.params[:v]
        else
            v_concat[1:dim1,1:dim1] = dist_in1.params[:v]
        end

        if dim2 == 1
            v_concat[end,end] = dist_in2.params[:v]
        else
            v_concat[dim1+1:end,dim1+1:end] = dist_in2.params[:v]
        end

        log_prior_pdf(x) = -0.5*(dim_tot*log(2pi) + log(det(v_concat)) + transpose(x-m_concat)*inv(v_concat)*(x-m_concat))

        function log_joint_dims(s::Array,dim1::Int64,dim2::Int64)
            if dim1 == 1
                if dim2 == 1
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1]::Number,s[end]::Number))
                else
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1],s[2:end]))
                end
            else
                if dim2 == 1
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1:dim1],s[end]))
                else
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1:dim1],s[dim1+1:end]))
                end
            end
        end

        log_joint(s) = log_joint_dims(s,dim1,dim2)

        #Optimization with gradient ascent
        d_log_joint(s) = ForwardDiff.gradient(log_joint, s)
        m_old = m_concat #initial point
        step_size = 0.01 #initial step size
        satisfied = 0
        step_count = 0
        m_total = zeros(dim_tot)
        m_average = zeros(dim_tot)
        m_new = zeros(dim_tot)
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
                if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < dim_tot*0.1
                    satisfied = 1
                end
            end
            if step_count > dim_tot*250
                satisfied = 1
            end
            m_old = m_new
        end
        m_post = m_new
        var_post = Hermitian(inv(- 1.0 .* ForwardDiff.jacobian(d_log_joint, m_post)))

        #decompose posterior estimations
        status[:count_update] = true

        if dim2 == 1
            mean2 = m_post[end]
            var2 = var_post[end]
            status[:message] = approxMessageBivariate(dist_in2.params[:m],dist_in2.params[:v],mean2,var2)
        else
            mean2 = m_post[dim1+1:end]
            var2 = var_post[dim1+1:end,dim1+1:end]
            status[:message] = approxMessageBivariate(dist_in2.params[:m],dist_in2.params[:v],mean2,var2)
        end

        if dim1 == 1
            mean1 = m_post[1]
            var1 = var_post[1]
            return approxMessageBivariate(dist_in1.params[:m],dist_in1.params[:v],mean1,var1)
        else
            mean1 = m_post[1:dim1]
            var1 = var_post[1:dim1,1:dim1]
            return approxMessageBivariate(dist_in1.params[:m],dist_in1.params[:v],mean1,var1)
        end
    end

end

function ruleSPBivariateLIn2MGN(msg_out::Message{Fout, Vout}, msg_in1::Message{F1, V1}, msg_in2::Message{F2, V2}, g::Function, status::Dict) where {Fout<:SoftFactor, Vout<:VariateType, F1<:Gaussian, V1<:VariateType, F2<:Gaussian, V2<:VariateType}

    if status[:updated]
        status[:updated] = false
        return status[:message]
    else
        dist_in1 = convert(ProbabilityDistribution{V1, GaussianMeanVariance}, msg_in1.dist)
        dist_in2 = convert(ProbabilityDistribution{V2, GaussianMeanVariance}, msg_in2.dist)

        m_concat = [dist_in1.params[:m];dist_in2.params[:m]]
        dim1 = dims(dist_in1)
        dim2 = dims(dist_in2)
        dim_tot = dim1 + dim2
        v_concat = zeros(dim_tot,dim_tot)

        if dim1 == 1
            v_concat[dim1,dim1] = dist_in1.params[:v]
        else
            v_concat[1:dim1,1:dim1] = dist_in1.params[:v]
        end

        if dim2 == 1
            v_concat[end,end] = dist_in2.params[:v]
        else
            v_concat[dim1+1:end,dim1+1:end] = dist_in2.params[:v]
        end

        log_prior_pdf(x) = -0.5*(dim_tot*log(2pi) + log(det(v_concat)) + transpose(x-m_concat)*inv(v_concat)*(x-m_concat))

        function log_joint_dims(s::Array,dim1::Int64,dim2::Int64)
            if dim1 == 1
                if dim2 == 1
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1]::Number,s[end]::Number))
                else
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1],s[2:end]))
                end
            else
                if dim2 == 1
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1:dim1],s[end]))
                else
                    return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1:dim1],s[dim1+1:end]))
                end
            end
        end

        log_joint(s) = log_joint_dims(s,dim1,dim2)
        #Optimization with gradient ascent
        d_log_joint(s) = ForwardDiff.gradient(log_joint, s)
        m_old = m_concat #initial point
        step_size = 0.01 #initial step size
        satisfied = 0
        step_count = 0
        m_total = zeros(dim_tot)
        m_average = zeros(dim_tot)
        m_new = zeros(dim_tot)
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
                if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < dim_tot*0.1
                    satisfied = 1
                end
            end
            if step_count > dim_tot*250
                satisfied = 1
            end
            m_old = m_new
        end
        m_post = m_new
        var_post = Hermitian(inv(- 1.0 .* ForwardDiff.jacobian(d_log_joint, m_post)))

        #decompose posterior estimations
        status[:updated] = true

        if dim1 == 1
            mean1 = m_post[1]
            var1 = var_post[1]
            status[:message] = approxMessageBivariate(dist_in1.params[:m],dist_in1.params[:v],mean1,var1)
        else
            mean1 = m_post[1:dim1]
            var1 = var_post[1:dim1,1:dim1]
            status[:message] = approxMessageBivariate(dist_in1.params[:m],dist_in1.params[:v],mean1,var1)
        end

        if dim2 == 1
            mean2 = m_post[end]
            var2 = var_post[end]
            return approxMessageBivariate(dist_in2.params[:m],dist_in2.params[:v],mean2,var2)
        else
            mean2 = m_post[dim1+1:end]
            var2 = var_post[dim1+1:end,dim1+1:end]
            return approxMessageBivariate(dist_in2.params[:m],dist_in2.params[:v],mean2,var2)
        end
    end

end

function ruleMBivariateLOutNGG(msg_out::Message{Fout, Vout}, msg_in1::Message{F1, V1}, msg_in2::Message{F2, V2}, g::Function, status::Dict) where {Fout<:SoftFactor, Vout<:VariateType, F1<:Gaussian, V1<:VariateType, F2<:Gaussian, V2<:VariateType}

    dist_in1 = convert(ProbabilityDistribution{V1, GaussianMeanVariance}, msg_in1.dist)
    dist_in2 = convert(ProbabilityDistribution{V2, GaussianMeanVariance}, msg_in2.dist)

    m_concat = [dist_in1.params[:m];dist_in2.params[:m]]

    dim1 = dims(dist_in1)
    dim2 = dims(dist_in2)
    dim_tot = dim1 + dim2
    v_concat = zeros(dim_tot,dim_tot)

    if dim1 == 1
        v_concat[dim1,dim1] = dist_in1.params[:v]
    else
        v_concat[1:dim1,1:dim1] = dist_in1.params[:v]
    end

    if dim2 == 1
        v_concat[end,end] = dist_in2.params[:v]
    else
        v_concat[dim1+1:end,dim1+1:end] = dist_in2.params[:v]
    end

    log_prior_pdf(x) = -0.5*(dim_tot*log(2pi) + log(det(v_concat)) + transpose(x-m_concat)*inv(v_concat)*(x-m_concat))

    function log_joint_dims(s::Array,dim1::Int64,dim2::Int64)
        if dim1 == 1
            if dim2 == 1
                return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1]::Number,s[end]::Number))
            else
                return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1],s[2:end]))
            end
        else
            if dim2 == 1
                return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1:dim1],s[end]))
            else
                return log_prior_pdf(s) + logPdf(msg_out.dist,g(s[1:dim1],s[dim1+1:end]))
            end
        end
    end

    log_joint(s) = log_joint_dims(s,dim1,dim2)

    #Optimization with gradient ascent
    d_log_joint(s) = ForwardDiff.gradient(log_joint, s)
    m_old = m_concat #initial point
    step_size = 0.01 #initial step size
    satisfied = 0
    step_count = 0
    m_total = zeros(dim_tot)
    m_average = zeros(dim_tot)
    m_new = zeros(dim_tot)
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
            if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < dim_tot*0.1
                satisfied = 1
            end
        end
        if step_count > dim_tot*250
            satisfied = 1
        end
        m_old = m_new
    end
    m_post = m_new
    var_post = inv(- 1.0 .* ForwardDiff.jacobian(d_log_joint, m_post))

    return ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_post, v=var_post)

end

#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::Bivariate{Laplace}, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    inbounds = Any[]
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
            haskey(interface_to_schedule_entry, node_interface) || error("This rule requires the incoming message on the out interface. Try altering execution order to ensure its availability.")
            if entry.message_update_rule == SPBivariateLOutNGG
                push!(inbounds, nothing)
            else
                push!(inbounds, interface_to_schedule_entry[inbound_interface])
            end
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
    status = "currentGraph().nodes[:$(node.id)].status"
    push!(inbounds, Dict{Symbol, Any}(:status => status,
                                      :keyword => false))

    return inbounds
end

#--------------------------
# Custom marginal inbounds collector
#--------------------------

function collectMarginalNodeInbounds(node::Bivariate, entry::MarginalEntry)
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

    # Push function and status to calling signature
    # The function needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))
    status = "currentGraph().nodes[:$(node.id)].status"
    push!(inbounds, Dict{Symbol, Any}(:status => status,
                                      :keyword => false))

    return inbounds
end
