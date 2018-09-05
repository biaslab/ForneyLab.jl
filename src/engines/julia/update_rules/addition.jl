export  
ruleSPAdditionOutVGG,
ruleSPAdditionOutVGP,
ruleSPAdditionOutVPG,
ruleSPAdditionOutVPP,
ruleSPAdditionIn1GVG,
ruleSPAdditionIn1PVG,
ruleSPAdditionIn1GVP,
ruleSPAdditionIn1PVP,
ruleSPAdditionIn2GGV,
ruleSPAdditionIn2PGV,
ruleSPAdditionIn2GPV,
ruleSPAdditionIn2PPV

function ruleSPAdditionOutVGG(
    msg_out::Nothing,
    msg_in1::Message{F1, V},
    msg_in2::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
    d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

    Message(V, GaussianMeanVariance, m=d_in1.params[:m] + d_in2.params[:m], v=d_in1.params[:v] + d_in2.params[:v])
end

function ruleSPAdditionIn2GGV(
    msg_out::Message{F1, V},
    msg_in1::Message{F2, V},
    msg_in2::Nothing) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
    d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

    Message(V, GaussianMeanVariance, m=d_out.params[:m] - d_in1.params[:m], v=d_out.params[:v] + d_in1.params[:v])
end

function ruleSPAdditionIn1GVG(
    msg_out::Message{F1, V},
    ::Nothing, 
    msg_in2::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn2GGV(msg_out, msg_in2, nothing)
end

function ruleSPAdditionOutVGP(  
    msg_out::Nothing,
    msg_in1::Message{F, V},
    msg_in2::Message{PointMass, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)

    Message(V, GaussianMeanVariance, m=d_in1.params[:m] + msg_in2.dist.params[:m], v=d_in1.params[:v])
end

function ruleSPAdditionOutVPG(
    ::Nothing, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Message{F, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionOutVGP(nothing, msg_in2, msg_in1)
end

function ruleSPAdditionIn1PVG(  
    msg_out::Message{PointMass, V},
    msg_in1::Nothing,
    msg_in2::Message{F, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

    Message(V, GaussianMeanVariance, m=msg_out.dist.params[:m] - d_in2.params[:m], v=d_in2.params[:v])
end

function ruleSPAdditionIn2PGV(
    msg_out::Message{PointMass, V}, 
    msg_in1::Message{F, V}, 
    msg_in2::Nothing) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}
    
    ruleSPAdditionIn1PVG(msg_out, nothing, msg_in1)
end

function ruleSPAdditionIn1GVP(  
    msg_out::Message{F, V},
    msg_in1::Nothing,
    msg_in2::Message{PointMass, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

    Message(V, GaussianMeanVariance, m=d_out.params[:m] - msg_in2.dist.params[:m], v=d_out.params[:v])
end

function ruleSPAdditionIn2GPV(
    msg_out::Message{F, V}, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Nothing) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

    ruleSPAdditionIn1GVP(msg_out, nothing, msg_in1)
end

function ruleSPAdditionOutVPP(
    msg_out::Nothing, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Message{PointMass, V}) where V<:Union{Univariate, Multivariate}

    Message(V, PointMass, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m])
end

function ruleSPAdditionIn2PPV(
    msg_out::Message{PointMass, V}, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Nothing) where V<:Union{Univariate, Multivariate}

    Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m])
end

function ruleSPAdditionIn1PVP(
    msg_out::Message{PointMass, V}, 
    msg_in1::Nothing, 
    msg_in2::Message{PointMass, V}) where V<:Union{Univariate, Multivariate}

    Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m])
end