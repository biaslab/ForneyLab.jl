export
ruleSPAdditionOutNGG,
ruleSPAdditionOutNGP,
ruleSPAdditionOutNPG,
ruleSPAdditionOutNPP,
ruleSPAdditionIn1GNG,
ruleSPAdditionIn1PNG,
ruleSPAdditionIn1GNP,
ruleSPAdditionIn1PNP,
ruleSPAdditionIn2GGN,
ruleSPAdditionIn2PGN,
ruleSPAdditionIn2GPN,
ruleSPAdditionIn2PPN,
ruleSPAdditionOutNSP,
ruleSPAdditionOutNPS,
ruleSPAdditionIn1PNS,
ruleSPAdditionIn1SNP,
ruleSPAdditionIn2PSN,
ruleSPAdditionIn2SPN,
ruleMAdditionNGG

function ruleSPAdditionOutNGG(
    msg_out::Nothing,
    msg_in1::Message{F1, V},
    msg_in2::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
    d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

    Message(V, GaussianMeanVariance, m=d_in1.params[:m] + d_in2.params[:m], v=d_in1.params[:v] + d_in2.params[:v])
end

function ruleSPAdditionIn2GGN(
    msg_out::Message{F1, V},
    msg_in1::Message{F2, V},
    msg_in2::Nothing) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
    d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

    Message(V, GaussianMeanVariance, m=d_out.params[:m] - d_in1.params[:m], v=d_out.params[:v] + d_in1.params[:v])
end

function ruleSPAdditionIn1GNG(
    msg_out::Message{F1, V},
    ::Nothing,
    msg_in2::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn2GGN(msg_out, msg_in2, nothing)
end

function ruleSPAdditionOutNGP(
    msg_out::Nothing,
    msg_in1::Message{F, V},
    msg_in2::Message{PointMass, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)

    Message(V, GaussianMeanVariance, m=d_in1.params[:m] + msg_in2.dist.params[:m], v=d_in1.params[:v])
end

function ruleSPAdditionOutNPG(
    ::Nothing,
    msg_in1::Message{PointMass, V},
    msg_in2::Message{F, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionOutNGP(nothing, msg_in2, msg_in1)
end

function ruleSPAdditionIn1PNG(
    msg_out::Message{PointMass, V},
    msg_in1::Nothing,
    msg_in2::Message{F, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

    Message(V, GaussianMeanVariance, m=msg_out.dist.params[:m] - d_in2.params[:m], v=d_in2.params[:v])
end

function ruleSPAdditionIn2PGN(
    msg_out::Message{PointMass, V},
    msg_in1::Message{F, V},
    msg_in2::Nothing) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn1PNG(msg_out, nothing, msg_in1)
end

function ruleSPAdditionIn1GNP(
    msg_out::Message{F, V},
    msg_in1::Nothing,
    msg_in2::Message{PointMass, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

    Message(V, GaussianMeanVariance, m=d_out.params[:m] - msg_in2.dist.params[:m], v=d_out.params[:v])
end

function ruleSPAdditionIn2GPN(
    msg_out::Message{F, V},
    msg_in1::Message{PointMass, V},
    msg_in2::Nothing) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn1GNP(msg_out, nothing, msg_in1)
end

function ruleSPAdditionOutNPP(
    msg_out::Nothing,
    msg_in1::Message{PointMass, V},
    msg_in2::Message{PointMass, V}) where V<:Union{Univariate, Multivariate}

    Message(V, PointMass, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m])
end

function ruleSPAdditionIn2PPN(
    msg_out::Message{PointMass, V},
    msg_in1::Message{PointMass, V},
    msg_in2::Nothing) where V<:Union{Univariate, Multivariate}

    Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m])
end

function ruleSPAdditionIn1PNP(
    msg_out::Message{PointMass, V},
    msg_in1::Nothing,
    msg_in2::Message{PointMass, V}) where V<:Union{Univariate, Multivariate}

    Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m])
end

function ruleSPAdditionOutNSP(
    msg_out::Nothing,
    msg_in1::Message{SampleList, Univariate},
    msg_in2::Message{PointMass, Univariate})

    new_s = msg_in2.dist.params[:m] .+ msg_in1.dist.params[:s]

    Message(Univariate, SampleList, s=new_s)
end

function ruleSPAdditionOutNSP(
    msg_out::Nothing,
    msg_in1::Message{SampleList, Multivariate},
    msg_in2::Message{PointMass, Multivariate})

    new_s = [msg_in2.dist.params[:m] .+ sample for sample in msg_in1.dist.params[:s]]

    Message(Multivariate, SampleList, s=new_s)
end

function ruleSPAdditionOutNPS(
    ::Nothing,
    msg_in1::Message{PointMass, V},
    msg_in2::Message{SampleList, V}) where {V<:Union{Univariate, Multivariate}}

    ruleSPAdditionOutNSP(nothing, msg_in2, msg_in1)
end

function ruleSPAdditionIn1PNS(
    msg_out::Message{PointMass, Univariate},
    msg_in1::Nothing,
    msg_in2::Message{SampleList, Univariate})

    new_s = msg_in2.dist.params[:m] .- msg_in1.dist.params[:s]

    Message(Univariate, SampleList, s=new_s)
end

function ruleSPAdditionIn1PNS(
    msg_out::Message{PointMass, Multivariate},
    msg_in1::Nothing,
    msg_in2::Message{SampleList, Multivariate})

    new_s = [msg_in2.dist.params[:m] .- sample for sample in msg_in1.dist.params[:s]]

    Message(Multivariate, SampleList, s=new_s)
end

function ruleSPAdditionIn2PSN(
    msg_out::Message{PointMass, V},
    msg_in1::Message{SampleList, V},
    msg_in2::Nothing) where {V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn1PNS(msg_out, nothing, msg_in1)
end

function ruleSPAdditionIn1SNP(
    msg_out::Message{SampleList, Univariate},
    msg_in1::Nothing,
    msg_in2::Message{PointMass, Univariate})

    new_s = -(msg_in2.dist.params[:m] .- msg_in1.dist.params[:s])

    Message(Univariate, SampleList, s=new_s)
end

function ruleSPAdditionIn1SNP(
    msg_out::Message{SampleList, Multivariate},
    msg_in1::Nothing,
    msg_in2::Message{PointMass, Multivariate})

    new_s = [-(msg_in2.dist.params[:m] .- sample) for sample in msg_in1.dist.params[:s]]

    Message(Multivariate, SampleList, s=new_s)
end

function ruleSPAdditionIn2SPN(
    msg_out::Message{SampleList, V},
    msg_in1::Message{PointMass, V},
    msg_in2::Nothing) where {V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn1SNP(msg_out, nothing, msg_in1)
end

function ruleMAdditionNGG(
    msg_out::Message{<:Gaussian, V}, 
    msg_in1::Message{<:Gaussian, V}, 
    msg_in2::Message{<:Gaussian, V}) where V<:Union{Univariate, Multivariate}

    d_out = convert(ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, msg_out.dist)
    d_in1 = convert(ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, msg_in1.dist)
    d_in2 = convert(ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, msg_in2.dist)

    xi_out = d_out.params[:xi]
    W_out = d_out.params[:w]
    xi_in1 = d_in1.params[:xi]
    W_in1 = d_in1.params[:w]
    xi_in2 = d_in2.params[:xi]
    W_in2 = d_in2.params[:w]

    return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[xi_in1+xi_out; xi_in2+xi_out], w=[W_in1+W_out W_out; W_out W_in2+W_out])
end
