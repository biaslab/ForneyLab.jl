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

function ruleSPAdditionOutVGG{F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}(
    msg_out::Void,
    msg_in1::Message{F1, V},
    msg_in2::Message{F2, V})

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
    d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

    Message(V, GaussianMeanVariance, m=d_in1.params[:m] + d_in2.params[:m], v=d_in1.params[:v] + d_in2.params[:v])
end

function ruleSPAdditionIn2GGV{F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}(
    msg_out::Message{F1, V},
    msg_in1::Message{F2, V},
    msg_in2::Void)

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
    d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

    Message(V, GaussianMeanVariance, m=d_out.params[:m] - d_in1.params[:m], v=d_out.params[:v] + d_in1.params[:v])
end

function ruleSPAdditionIn1GVG{F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}(
    msg_out::Message{F1, V},
    ::Void, 
    msg_in2::Message{F2, V})

    ruleSPAdditionIn2GGV(msg_out, msg_in2, nothing)
end

function ruleSPAdditionOutVGP{F<:Gaussian, V<:Union{Univariate, Multivariate}}(  
    msg_out::Void,
    msg_in1::Message{F, V},
    msg_in2::Message{PointMass, V})

    d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)

    Message(V, GaussianMeanVariance, m=d_in1.params[:m] + msg_in2.dist.params[:m], v=d_in1.params[:v])
end

function ruleSPAdditionOutVPG{F<:Gaussian, V<:Union{Univariate, Multivariate}}(
    ::Void, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Message{F, V})

    ruleSPAdditionOutVGP(nothing, msg_in2, msg_in1)
end

function ruleSPAdditionIn1PVG{F<:Gaussian, V<:Union{Univariate, Multivariate}}(  
    msg_out::Message{PointMass, V},
    msg_in1::Void,
    msg_in2::Message{F, V})

    d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

    Message(V, GaussianMeanVariance, m=msg_out.dist.params[:m] - d_in2.params[:m], v=d_in2.params[:v])
end

function ruleSPAdditionIn2PGV{F<:Gaussian, V<:Union{Univariate, Multivariate}}(
    msg_out::Message{PointMass, V}, 
    msg_in1::Message{F, V}, 
    msg_in2::Void)
    
    ruleSPAdditionIn1PVG(msg_out, nothing, msg_in1)
end

function ruleSPAdditionIn1GVP{F<:Gaussian, V<:Union{Univariate, Multivariate}}(  
    msg_out::Message{F, V},
    msg_in1::Void,
    msg_in2::Message{PointMass, V})

    d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

    Message(V, GaussianMeanVariance, m=d_out.params[:m] - msg_in2.dist.params[:m], v=d_out.params[:v])
end

function ruleSPAdditionIn2GPV{F<:Gaussian, V<:Union{Univariate, Multivariate}}(
    msg_out::Message{F, V}, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Void)

    ruleSPAdditionIn1GVP(msg_out, nothing, msg_in1)
end

function ruleSPAdditionOutVPP{V<:Union{Univariate, Multivariate}}(
    msg_out::Void, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Message{PointMass, V})

    Message(V, PointMass, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m])
end

function ruleSPAdditionIn2PPV{V<:Union{Univariate, Multivariate}}(
    msg_out::Message{PointMass, V}, 
    msg_in1::Message{PointMass, V}, 
    msg_in2::Void)

    Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m])
end

function ruleSPAdditionIn1PVP{V<:Union{Univariate, Multivariate}}(
    msg_out::Message{PointMass, V}, 
    msg_in1::Void, 
    msg_in2::Message{PointMass, V})

    Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m])
end