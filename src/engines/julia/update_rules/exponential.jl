export
ruleSPExponentialOutNG,
ruleSPExponentialOutNP,
ruleSPExponentialIn1LN,
ruleSPExponentialIn1PN

function ruleSPExponentialOutNG(msg_out::Nothing, 
								msg_in1::Message{F, Univariate}) where F<:Gaussian

    d_in1 = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_in1.dist)

    return Message(Univariate, LogNormal, m=d_in1.params[:m], s=d_in1.params[:v])
end

ruleSPExponentialIn1LN(msg_out::Message{LogNormal, Univariate}, msg_in1::Nothing) = 
	Message(Univariate, GaussianMeanVariance, m=msg_out.dist.params[:m], v=msg_out.dist.params[:s])

ruleSPExponentialOutNP(msg_out::Nothing, msg_in1::Message{PointMass, Univariate}) = 
	Message(Univariate, PointMass, m=exp(msg_in1.dist.params[:m]))

ruleSPExponentialIn1PN(msg_out::Message{PointMass, Univariate}, msg_in1::Nothing) = 
	Message(Univariate, PointMass, m=log(msg_out.dist.params[:m]))