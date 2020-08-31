export
ruleSPBetaOutNMM,
ruleSPBetaMNM,
ruleSPBetaMMN,
ruleVBBetaOut

ruleSPBetaOutNMM(msg_out::Nothing, msg_a::Message{PointMass, Univariate}, msg_b::Message{PointMass, Univariate}) = Message(Univariate, Beta, a=deepcopy(msg_a.dist.params[:m]), b=deepcopy(msg_b.dist.params[:m]))

function ruleSPBetaOutNMM(msg_out::Nothing, msg_a::Message{F1, Univariate}, msg_b::Message{F2, Univariate}) where {F1<:Union{Function, FactorNode},F2<:Union{Function, FactorNode}}
    n_samples = default_n_samples
    samples_a = sample(msg_a.dist,n_samples)
    samples_b = sample(msg_b.dist,n_samples)
    s_list = Vector{Float64}(undef, n_samples)
    for n=1:n_samples
        c = ProbabilityDistribution(Univariate,Beta,a=samples_a[n],b=samples_b[n])
        s_list[n] = sample(c)
    end
    w_list = ones(n_samples)/n_samples
    Message(Univariate,SampleList,s=s_list,w=w_list)
end

#will not work for incoming function message for now. we need to write a sampler for function messages
#1000 samples is chosen arbitrarily. Better to be defined by the user.
function ruleSPBetaMNM(msg_out::Message{F1,Univariate}, msg_a::Nothing, msg_b::Message{F2,Univariate}) where {F1<:Union{Function, FactorNode},F2<:Union{Function, FactorNode}}
    n_samples = default_n_samples
    samples_out = sample(msg_out.dist,n_samples)
    samples_b = sample(msg_b.dist,n_samples)
    logp(x,a,b) = (a-1)*log(x) + (b-1)*log(1.0-x) - loggamma(a) - loggamma(b) + loggamma(a+b)
    logp(a) = sum(logp.(samples_out,a,samples_b))/n_samples
    Message(Univariate, Function, log_pdf = logp)
end

function ruleSPBetaMMN(msg_out::Message{F1,Univariate}, msg_a::Message{F2,Univariate}, msg_b::Nothing) where {F1<:Union{Function, FactorNode},F2<:Union{Function, FactorNode}}
    n_samples = default_n_samples
    samples_out = sample(msg_out.dist,n_samples)
    samples_a = sample(msg_a.dist,n_samples)
    logp(x,a,b) = (a-1)*log(x) + (b-1)*log(1.0-x) - loggamma(a) - loggamma(b) + loggamma(a+b)
    logp(b) = sum(logp.(samples_out,samples_a,b))/n_samples
    Message(Univariate, Function, log_pdf = logp)
end

ruleVBBetaOut(marg_out::Any, dist_a::ProbabilityDistribution{Univariate}, dist_b::ProbabilityDistribution{Univariate, PointMass}) = Message(Univariate, Beta, a=unsafeMean(dist_a), b=unsafeMean(dist_b))
