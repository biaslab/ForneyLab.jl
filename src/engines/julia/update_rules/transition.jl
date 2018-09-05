export
ruleSPTransitionOutVCP,
ruleSPTransitionIn1CVP,
ruleVBTransitionOut,
ruleVBTransitionIn1,
ruleVBTransitionA,
ruleSVBTransitionOutVCD,
ruleSVBTransitionIn1CVD,
ruleSVBTransitionADV,
ruleMTransitionCCD

function ruleSPTransitionOutVCP(msg_out::Nothing,
                                msg_in1::Message{Categorical, Univariate},
                                msg_a::Message{PointMass, MatrixVariate})

    a = msg_a.dist.params[:m]*msg_in1.dist.params[:p]
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionIn1CVP(msg_out::Message{Categorical, Univariate},
                                msg_in1::Nothing,
                                msg_a::Message{PointMass, MatrixVariate})

    a = msg_a.dist.params[:m]'*msg_out.dist.params[:p]
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionOut(   dist_out::Nothing,
                                dist_in1::ProbabilityDistribution,
                                dist_a::ProbabilityDistribution{MatrixVariate})

    a = exp.(unsafeLogMean(dist_a)*unsafeMeanVector(dist_in1))
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionIn1(   dist_out::ProbabilityDistribution,
                                dist_in1::Nothing,
                                dist_a::ProbabilityDistribution{MatrixVariate})

    a = exp.(unsafeLogMean(dist_a)'*unsafeMeanVector(dist_out))
    
    Message(Univariate, Categorical, p=a./sum(a))
end

ruleVBTransitionA(  dist_out::ProbabilityDistribution,
                    dist_in1::ProbabilityDistribution,
                    dist_a::Nothing) =
    Message(MatrixVariate, Dirichlet, a=unsafeMeanVector(dist_out)*unsafeMeanVector(dist_in1)' .+ 1)

function ruleSVBTransitionOutVCD(   dist_out::Nothing,
                                    msg_in1::Message{Categorical, Univariate},
                                    dist_a::ProbabilityDistribution{MatrixVariate})

    a = exp.(unsafeLogMean(dist_a))*msg_in1.dist.params[:p]

    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSVBTransitionIn1CVD(   msg_out::Message{Categorical, Univariate},
                                    dist_in1::Nothing,
                                    dist_a::ProbabilityDistribution{MatrixVariate})

    a = exp.(unsafeLogMean(dist_a))'*msg_out.dist.params[:p]

    Message(Univariate, Categorical, p=a./sum(a))
end

ruleSVBTransitionADV(   dist_out_in1::ProbabilityDistribution{Multivariate, Contingency},
                        dist_a::Nothing)=
    Message(MatrixVariate, Dirichlet, a=dist_out_in1.params[:p] .+ 1)

function ruleMTransitionCCD(msg_out::Message{Categorical, Univariate},
                            msg_in1::Message{Categorical, Univariate},
                            dist_a::ProbabilityDistribution{MatrixVariate})

    B = Diagonal(msg_out.dist.params[:p])*exp.(unsafeLogMean(dist_a))*Diagonal(msg_in1.dist.params[:p])

    ProbabilityDistribution(Multivariate, Contingency, p=B./sum(B))
end
