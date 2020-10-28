export
ruleSPTransitionOutNPP,
ruleSPTransitionIn1PNP,
ruleSPTransitionOutNCP,
ruleSPTransitionIn1CNP,
ruleVBTransitionOut,
ruleVBTransitionIn1,
ruleVBTransitionA,
ruleSVBTransitionOutVCD,
ruleSVBTransitionIn1CVD,
ruleSVBTransitionADV,
ruleMTransitionCCD,
ruleMTransitionCCN

# Note that vanishing probabilities are softened to prevent singularities
function ruleSPTransitionOutNPP(msg_out::Nothing,
                                msg_in1::Message{PointMass, Multivariate},
                                msg_a::Message{PointMass, MatrixVariate})

    a = clamp.(msg_a.dist.params[:m]*msg_in1.dist.params[:m], tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionIn1PNP(msg_out::Message{PointMass, Multivariate},
                                msg_in1::Nothing,
                                msg_a::Message{PointMass, MatrixVariate})

    a = clamp.(msg_a.dist.params[:m]'*msg_out.dist.params[:m], tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionOutNCP(msg_out::Nothing,
                                msg_in1::Message{Categorical, Univariate},
                                msg_a::Message{PointMass, MatrixVariate})

    a = clamp.(msg_a.dist.params[:m]*msg_in1.dist.params[:p], tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionIn1CNP(msg_out::Message{Categorical, Univariate},
                                msg_in1::Nothing,
                                msg_a::Message{PointMass, MatrixVariate})

    a = clamp.(msg_a.dist.params[:m]'*msg_out.dist.params[:p], tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionOut(   dist_out::Nothing,
                                dist_in1::ProbabilityDistribution,
                                dist_a::ProbabilityDistribution{MatrixVariate})

    a = clamp.(exp.(unsafeLogMean(dist_a)*unsafeMeanVector(dist_in1)), tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionIn1(   dist_out::ProbabilityDistribution,
                                dist_in1::Nothing,
                                dist_a::ProbabilityDistribution{MatrixVariate})

    a = clamp.(exp.(unsafeLogMean(dist_a)'*unsafeMeanVector(dist_out)), tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

ruleVBTransitionA(  dist_out::ProbabilityDistribution,
                    dist_in1::ProbabilityDistribution,
                    dist_a::Nothing) =
    Message(MatrixVariate, Dirichlet, a=unsafeMeanVector(dist_out)*unsafeMeanVector(dist_in1)' .+ 1)

function ruleSVBTransitionOutVCD(   dist_out::Nothing,
                                    msg_in1::Message{Categorical, Univariate},
                                    dist_a::ProbabilityDistribution{MatrixVariate})

    a = clamp.(exp.(unsafeLogMean(dist_a))*msg_in1.dist.params[:p], tiny, Inf)

    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSVBTransitionIn1CVD(   msg_out::Message{Categorical, Univariate},
                                    dist_in1::Nothing,
                                    dist_a::ProbabilityDistribution{MatrixVariate})

    a = clamp.(exp.(unsafeLogMean(dist_a))'*msg_out.dist.params[:p], tiny, Inf)

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

function ruleMTransitionCCN(msg_out::Message{Categorical, Univariate},
                            msg_in1::Message{Categorical, Univariate},
                            msg_a::Message{PointMass, MatrixVariate})

    B = Diagonal(msg_out.dist.params[:p])*msg_a.dist.params[:m]*Diagonal(msg_in1.dist.params[:p])

    ProbabilityDistribution(Multivariate, Contingency, p=B./sum(B))
end
