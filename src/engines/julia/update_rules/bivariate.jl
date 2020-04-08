export
ruleSPBivariateLIn1MNG,
ruleSPBivariateLIn2MGN,
ruleSPBivariateLOutNGG,
ruleMBivariateLOutNGG

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F, Univariate}, msg_in2::Message{F, Univariate}, g::Function, status::Dict) where {F<:Gaussian}
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

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F, Multivariate}, msg_in2::Message{F, Univariate}, g::Function, status::Dict) where {F<:Gaussian}
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

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F, Univariate}, msg_in2::Message{F, Multivariate}, g::Function, status::Dict) where {F<:Gaussian}
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

function ruleSPBivariateLOutNGG(msg_out::Nothing, msg_in1::Message{F, Multivariate}, msg_in2::Message{F, Multivariate}, g::Function, status::Dict) where {F<:Gaussian}
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

KL_between_2Gaussians(m1::Number, v1::Number, m2::Number, v2::Number) = log(sqrt(v2)/sqrt(v1)) + (v1+(m1-m2)^2)/(2*v2) - 0.5

function KL_between_2Gaussians(m1::Array, v1, m2::Array, v2)
    v2_inv = inv(v2+1e-5*diageye(length(m1)))
    return 0.5*(log(det(v2)/det(v1)) - length(m1) + tr(v2_inv*v1) + transpose(m2-m1)*v2_inv*(m2-m1))
end

function approxMessageBivariate(m_prior::Number,v_prior::Number,m_post::Number,v_post::Number)
    # #gradient descent to find m_message, v_message that minimizes KL objective
    # m_prod(m1,v1,m2,v2) = (v1*m2 + v2*m1)/(v1+v2) #mean of prod of two univariate Gaussian dist.s
    # v_prod(m1,v1,m2,v2) = 1/(1/v1+1/v2) #var of prod of two univariate Gaussian dist.s
    #
    # KL_objective(m,v) = 0.5*(KL_between_2Gaussians(m_prod(m,v,m_prior,v_prior),v_prod(m,v,m_prior,v_prior),m_post,v_post)+KL_between_2Gaussians(m_post,v_post,m_prod(m,v,m_prior,v_prior),v_prod(m,v,m_prior,v_prior)))
    # #KL_objective(m,v) = KL_between_2Gaussians(m_prod(m,v,m_prior,v_prior),v_prod(m,v,m_prior,v_prior),m_post,v_post)
    #
    # step_size = 0.01
    # satisfied = 0
    # step_count = 0
    # m_message, s_message =  m_post, sqrt(v_post) #initial values for mean and std
    # m_total, s_total = 0, 0
    # m_average, s_average = 0, 0
    # m_new, s_new = 0, 0
    # @show m_post
    # @show v_post
    # while satisfied == 0
    #     KL_m(m) = KL_objective(m,s_message^2)
    #     KL_s(s) = KL_objective(m_message,s^2)
    #     d_KL_m(m) = ForwardDiff.derivative(KL_m,m)
    #     d_KL_s(s) = ForwardDiff.derivative(KL_s,s)
    #     m_new = m_message - step_size*d_KL_m(m_message)
    #     @show d_KL_s(s_message)
    #     s_new = s_message - step_size*d_KL_s(s_message)
    #     if KL_objective(m_new,s_new^2) < KL_objective(m_message,s_message^2)
    #         proposal_step_size = 10*step_size
    #         m_proposal = m_message - proposal_step_size*d_KL_m(m_message)
    #         s_proposal = s_message - proposal_step_size*d_KL_s(s_message)
    #         if KL_objective(m_proposal,s_proposal^2) < KL_objective(m_new,s_new^2)
    #             m_new = m_proposal
    #             s_new = s_proposal
    #             step_size = proposal_step_size
    #         end
    #     else
    #         step_size = 0.1*step_size
    #         m_new = m_message - step_size*d_KL_m(m_message)
    #         s_new = s_message - step_size*d_KL_s(s_message)
    #     end
    #     step_count += 1
    #     m_total += m_message
    #     m_average = m_total / step_count
    #     s_total += s_message
    #     s_average = s_total / step_count
    #     if step_count > 10
    #         if abs((m_new-m_average)/m_average) < 0.001 && abs((s_new-s_average)/s_average) < 0.001
    #             satisfied = 1
    #             @show m_message
    #             @show s_message^2
    #         end
    #     end
    #     if step_count > 250
    #         satisfied = 1
    #         @show m_message
    #         @show s_message^2
    #     end
    #     m_message = m_new
    #     s_message = s_new
    # end
    # v_message = s_message^2
    # @show m_message
    # @show s_message^2
    # return Message(Univariate, GaussianMeanVariance, m=m_message, v=v_message)
    if abs(v_prior-v_post) < 1e-5
        v_message = 1e-5
    else
        v_message = v_prior*v_post/(v_prior-v_post)
    end
    m_message = (m_post*(v_prior+v_message) - m_prior*v_message)/v_prior
    return Message(Univariate, GaussianMeanVariance, m=m_message, v=v_message)
end

function approxMessageBivariate(m_prior::Array,v_prior::Array,m_post::Array,v_post::Array)
    # #gradient descent to find m_message, v_message that minimizes KL objective
    # #Mean and Variance of prod of two multivariate Gaussian dist.s
    # function G_prod(m1,v1,m2,v2)
    #     L = cholesky(Hermitian(v1+v2)).L
    #     inv_L = inv(L)
    #     V1, V2 = inv_L*v1, inv_L*v2
    #     M1, M2 = inv_L*m1, inv_L*m2
    #     return transpose(V2)*m1 + transpose(V1)*m2, Hermitian(transpose(V1)*V2)
    # end
    #
    # function KL_objective(m,v)
    #     m_p, v_p = G_prod(m,v,m_prior,v_prior)
    #     0.5*(KL_between_2Gaussians(m_p,v_p,m_post,v_post)+KL_between_2Gaussians(m_post,v_post,m_p,v_p))
    # end
    #
    # step_size = 0.001
    # satisfied = 0
    # step_count = 0
    # m_message, s_message =  m_post, cholesky(v_post).L #initial values for mean and lower triangular cholesky component
    # m_total, s_total = zeros(length(m_post)), zeros(length(m_post),length(m_post))
    # m_average, s_average = zeros(length(m_post)), zeros(length(m_post),length(m_post))
    # m_new, s_new = zeros(length(m_post)), zeros(length(m_post),length(m_post))
    # while satisfied == 0
    #     KL_m(m) = KL_objective(m,s_message*transpose(s_message))
    #     KL_s(s) = KL_objective(m_message,s*transpose(s))
    #     d_KL_m(m) = ForwardDiff.gradient(KL_m,m)
    #     d_KL_s(s) = ForwardDiff.gradient(KL_s,s)
    #     #be sure that v_new is positive definite
    #     p_satisfied = 0
    #     while p_satisfied == 0
    #         m_new = m_message .- step_size.*d_KL_m(m_message)
    #         s_new = s_message .- step_size.*d_KL_s(s_message)
    #         if all(diag(s_new).>0)
    #             p_satisfied = 1
    #         else
    #             step_size = 0.1*step_size
    #         end
    #     end
    #     if KL_objective(m_new,s_new*transpose(s_new)) < KL_objective(m_message,s_message*transpose(s_message))
    #         proposal_step_size = 10*step_size
    #         m_proposal = m_message .- proposal_step_size.*d_KL_m(m_message)
    #         s_proposal = s_message .- proposal_step_size.*d_KL_s(s_message)
    #         if all(diag(s_proposal).>0) #be sure that v_proposal is positive definite
    #             if KL_objective(m_proposal,s_proposal*transpose(s_proposal)) < KL_objective(m_new,s_new*transpose(s_new))
    #                 m_new = m_proposal
    #                 s_new = s_proposal
    #                 step_size = proposal_step_size
    #             end
    #         end
    #     else
    #         step_size = 0.1*step_size
    #         m_new = m_message .- step_size.*d_KL_m(m_message)
    #         s_new = s_message .- step_size.*d_KL_s(s_message)
    #     end
    #     step_count += 1
    #     m_total .+= m_message
    #     m_average = m_total ./ step_count
    #     s_total .+= s_message
    #     s_average = s_total ./ step_count
    #     if step_count > 10
    #         if sum(sqrt.(((m_new.-m_average)./m_average).^2)) < length(m_post)*0.1 && sum(sqrt.(((s_new.-s_average)./s_average).^2)) < length(m_post)*0.1
    #             satisfied = 1
    #         end
    #     end
    #     if step_count > length(m_post)*250
    #         satisfied = 1
    #     end
    #     m_message = m_new
    #     s_message = s_new
    #     @show s_message
    # end
    # v_message = s_message*transpose(s_message)
    # return Message(Multivariate, GaussianMeanVariance, m=m_message, v=v_message)
    w_prior, w_post = inv(v_prior+2e-5*diageye(length(m_prior))), inv(v_post+1e-5*diageye(length(m_prior)))
    w_message = w_post - w_prior
    xi_message = (w_prior+w_message)*m_post - w_prior*m_prior
    return Message(Multivariate, GaussianWeightedMeanPrecision, xi=xi_message, w=w_message)
end


function ruleSPBivariateLIn1MNG(msg_out::Message{Fout, Vout}, msg_in1::Message{F1, V1}, msg_in2::Message{F2, V2}, g::Function, status::Dict) where {Fout<:SoftFactor, Vout<:VariateType, F1<:Gaussian, V1<:VariateType, F2<:Gaussian, V2<:VariateType}

    if status[:count_update] == 1
        status[:count_update] = 0
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
        status[:count_update] = 1

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

function ruleSPBivariateLIn2MGN(msg_out::Message{Fout, Vout}, msg_in1::Message{F1, V1}, msg_in2::Message{F1, V2}, g::Function, status::Dict) where {Fout<:SoftFactor, Vout<:VariateType, F1<:Gaussian, V1<:VariateType, F2<:Gaussian, V2<:VariateType}

    if status[:count_update] == 1
        status[:count_update] = 0
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
        status[:count_update] = 1

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

function ruleMBivariateLOutNGG(msg_out::Message{Fout, Vout}, msg_in1::Message{F, V1}, msg_in2::Message{F, V2}, g::Function, status::Dict) where {Fout<:SoftFactor, Vout<:VariateType, F<:Gaussian, V1<:VariateType, V2<:VariateType}

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
