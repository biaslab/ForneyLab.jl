export
ruleSPBetaOutNMM,
ruleSPBetaMNM,
ruleSPBetaMMN,
ruleVBBetaOut,
ruleVBBetaA,
ruleVBBetaB

ruleSPBetaOutNMM(msg_out::Nothing, msg_a::Message{PointMass, Univariate}, msg_b::Message{PointMass, Univariate}) = Message(Univariate, Beta, a=deepcopy(msg_a.dist.params[:m]), b=deepcopy(msg_b.dist.params[:m]))

function ruleSPBetaOutNMM(msg_out::Nothing, msg_a::Message{F1, Univariate}, msg_b::Message{F2, Univariate}) where {F1<:FactorNode, F2<:FactorNode}
    n_samples = default_n_samples
    samples_a = sample(msg_a.dist,n_samples)
    samples_b = sample(msg_b.dist,n_samples)
    s_list = betainvcdf.(samples_a, samples_b, rand(n_samples))
    w_list = ones(n_samples)/n_samples
    Message(Univariate,SampleList,s=s_list,w=w_list)
end

# 1000 (default_n_samples) samples is chosen arbitrarily. Better to be defined by the user.
function ruleSPBetaMNM(msg_out::Message{F1,Univariate}, msg_a::Nothing, msg_b::Message{F2,Univariate}) where {F1<:FactorNode, F2<:FactorNode}
    n_samples = default_n_samples
    samples_out = sample(msg_out.dist,n_samples)
    samples_b = sample(msg_b.dist,n_samples)
    logp(x,a,b) = (a-1)*log(x) + (b-1)*log(1.0-x) - loggamma(a) - loggamma(b) + loggamma(a+b)
    logp(a) = sum(logp.(samples_out,a,samples_b))/n_samples
    Message(Univariate, Function, log_pdf = logp)
end

function ruleSPBetaMMN(msg_out::Message{F1,Univariate}, msg_a::Message{F2,Univariate}, msg_b::Nothing) where {F1<:FactorNode, F2<:FactorNode}
    n_samples = default_n_samples
    samples_out = sample(msg_out.dist,n_samples)
    samples_a = sample(msg_a.dist,n_samples)
    logp(x,a,b) = (a-1)*log(x) + (b-1)*log(1.0-x) - loggamma(a) - loggamma(b) + loggamma(a+b)
    logp(b) = sum(logp.(samples_out,samples_a,b))/n_samples
    Message(Univariate, Function, log_pdf = logp)
end

ruleVBBetaOut(marg_out::Any, dist_a::ProbabilityDistribution{Univariate}, dist_b::ProbabilityDistribution{Univariate}) = Message(Univariate, Beta, a=unsafeMean(dist_a), b=unsafeMean(dist_b))

function ruleVBBetaA( dist_out::ProbabilityDistribution{Univariate},
    dist_a::Any, dist_b::ProbabilityDistribution{Univariate})

    logGamAsumB(a) = sum(loggamma.(a.+sample(dist_b,default_n_samples)))/default_n_samples
    logp(a) = (a-1)*unsafeLogMean(dist_out) + logGamAsumB(a) - loggamma(a)

    Message(Univariate, Function, log_pdf = logp)
end

function ruleVBBetaB( dist_out::ProbabilityDistribution{Univariate},
    dist_a::ProbabilityDistribution{Univariate}, dist_b::Any)

    logGamAsumB(b) = sum(loggamma.(b.+sample(dist_a,default_n_samples)))/default_n_samples
    logp(b) = (b-1)*unsafeMirroredLogMean(dist_out) + logGamAsumB(b) - loggamma(b)

    Message(Univariate, Function, log_pdf = logp)
end
