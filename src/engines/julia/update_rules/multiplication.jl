export
ruleSPMultiplicationOutNPG,
ruleSPMultiplicationOutNGP,
ruleSPMultiplicationOutNPP,
ruleSPMultiplicationIn1GNP,
ruleSPMultiplicationIn1PNP,
ruleSPMultiplicationAGPN,
ruleSPMultiplicationAPPN

#-------------------------------
# Univariate (in1 and a commute)
#-------------------------------

function ruleSPMultiplicationOutNPG(msg_out::Nothing,
                                    msg_in1::Message{PointMass, Univariate},
                                    msg_a::Message{F, Univariate}) where F<:Gaussian

    d_a = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_a.dist)

    Message(Univariate, GaussianMeanVariance, m=msg_in1.dist.params[:m]*d_a.params[:m], v=d_a.params[:v]*msg_in1.dist.params[:m]^2)
end

ruleSPMultiplicationOutNGP(msg_out::Nothing, msg_in1::Message{F, Univariate}, msg_a::Message{PointMass, Univariate}) where F<:Gaussian = 
    ruleSPMultiplicationOutNPG(nothing, msg_a, msg_in1)

function ruleSPMultiplicationIn1GNP(msg_out::Message{F, Univariate},
                                    msg_in1::Nothing,
                                    msg_a::Message{PointMass, Univariate}) where F<:Gaussian

    d_out = convert(ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}, msg_out.dist)
    
    a = msg_a.dist.params[:m]

    Message(Univariate, GaussianWeightedMeanPrecision, xi=a*d_out.params[:xi], w=d_out.params[:w]*a^2)
end

ruleSPMultiplicationAGPN(msg_out::Message{F, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_a::Nothing) where F<:Gaussian = 
    ruleSPMultiplicationIn1GNP(msg_out, nothing, msg_in1)

ruleSPMultiplicationOutNPP(msg_out::Nothing, msg_in1::Message{PointMass, Univariate}, msg_a::Message{PointMass, Univariate}) = 
    Message(Univariate, PointMass, m=msg_in1.dist.params[:m]*msg_a.dist.params[:m])

ruleSPMultiplicationIn1PNP(msg_out::Message{PointMass, Univariate}, msg_in1::Nothing, msg_a::Message{PointMass, Univariate}) = 
    Message(Univariate, PointMass, m=msg_out.dist.params[:m]/msg_a.dist.params[:m])

ruleSPMultiplicationAPPN(msg_out::Message{PointMass, Univariate}, msg_in1::Message{PointMass, Univariate}, msg_a::Nothing) = 
    Message(Univariate, PointMass, m=msg_out.dist.params[:m]/msg_in1.dist.params[:m])


#--------------------------------------------
# Univariate*multivariate (in1 and a commute)
#--------------------------------------------

function ruleSPMultiplicationOutNPG(msg_out::Nothing,
                                    msg_in1::Message{PointMass, Multivariate},
                                    msg_a::Message{F, Univariate}) where F<:Gaussian

    d_a = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_a.dist)

    Message(Multivariate, GaussianMeanVariance, m=msg_in1.dist.params[:m]*d_a.params[:m], v=d_a.params[:v]*msg_in1.dist.params[:m]*msg_in1.dist.params[:m]')
end

ruleSPMultiplicationOutNGP(msg_out::Nothing, msg_in1::Message{F, Univariate}, msg_a::Message{PointMass, Multivariate}) where F<:Gaussian = 
    ruleSPMultiplicationOutNPG(nothing, msg_a, msg_in1)

function ruleSPMultiplicationOutNGP(msg_out::Nothing, 
                                    msg_in1::Message{F, Multivariate}, 
                                    msg_a::Message{PointMass, Univariate}) where F<:Gaussian

    d_in1 = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in1.dist)

    Message(Multivariate, GaussianMeanVariance, m=d_in1.params[:m]*msg_a.dist.params[:m], v=d_in1.params[:v]*msg_a.dist.params[:m]^2)
end

ruleSPMultiplicationOutNPG(msg_out::Nothing, msg_in1::Message{PointMass, Univariate}, msg_a::Message{F, Multivariate}) where F<:Gaussian = 
    ruleSPMultiplicationOutNGP(nothing, msg_a, msg_in1)


#------------------------------------------------------
# MatrixVariate*multivariate (in1 and a do NOT commute)
#------------------------------------------------------

function ruleSPMultiplicationOutNGP(msg_out::Nothing,
                                    msg_in1::Message{F, Multivariate},
                                    msg_a::Message{PointMass, MatrixVariate}) where F<:Gaussian

    d_in1 = convert(ProbabilityDistribution{Multivariate, GaussianMeanVariance}, msg_in1.dist)

    A = msg_a.dist.params[:m]

    Message(Multivariate, GaussianMeanVariance, m=A*d_in1.params[:m], v=A*d_in1.params[:v]*A')
end

function ruleSPMultiplicationIn1GNP(msg_out::Message{F, Multivariate},
                                    msg_in1::Nothing,
                                    msg_a::Message{PointMass, MatrixVariate}) where F<:Gaussian

    d_out = convert(ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}, msg_out.dist)

    A = msg_a.dist.params[:m]
    W = A'*d_out.params[:w]*A
    W = W + tiny*diageye(size(W)[1]) # Ensure precision is always invertible

    Message(Multivariate, GaussianWeightedMeanPrecision, xi=A'*d_out.params[:xi], w=W)
end

ruleSPMultiplicationOutNPP(msg_out::Nothing, msg_in1::Message{PointMass, Multivariate}, msg_a::Message{PointMass, MatrixVariate}) = 
    Message(Multivariate, PointMass, m=msg_a.dist.params[:m]*msg_in1.dist.params[:m])

ruleSPMultiplicationIn1PNP(msg_out::Message{PointMass, Multivariate}, msg_in1::Nothing, msg_a::Message{PointMass, MatrixVariate}) = 
    Message(Multivariate, PointMass, m=pinv(msg_a.dist.params[:m])*msg_out.dist.params[:m])