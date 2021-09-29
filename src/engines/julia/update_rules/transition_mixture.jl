export
ruleSPTransitionMixtureOutNCCPX,
ruleSPTransitionMixtureIn1CNCPX,
ruleSPTransitionMixtureZCCNPX,
ruleSVBTransitionMixtureOutNCCDX,
ruleSVBTransitionMixtureIn1CNCDX,
ruleSVBTransitionMixtureZCCNDX, 
ruleSVBTransitionMixtureA,
ruleMTransitionMixtureCCCDX,
ruleMTransitionMixtureCCCNX

function ruleSPTransitionMixtureOutNCCPX(msg_out::Nothing,
                                         msg_in1::Message{Categorical},
                                         msg_switch::Message{Categorical},
                                         msg_factors::Vararg{Message{PointMass, MatrixVariate}})
    z_bar = msg_switch.dist.params[:p]
    n_factors = length(msg_factors)
    d = dims(msg_factors[1].dist)[1] # Dimensionality of output
    a = tiny*ones(d)
    for k=1:n_factors
        a += z_bar[k]*msg_factors[k].dist.params[:m]*msg_in1.dist.params[:p]
    end
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionMixtureIn1CNCPX(msg_out::Message{Categorical},
                                         msg_in1::Nothing,
                                         msg_switch::Message{Categorical},
                                         msg_factors::Vararg{Message{PointMass, MatrixVariate}})
    z_bar = msg_switch.dist.params[:p]
    n_factors = length(msg_factors)
    d = dims(msg_factors[1].dist)[2] # Dimensionality of input
    a = tiny*ones(d)
    for k=1:n_factors
        a += z_bar[k]*msg_factors[k].dist.params[:m]'*msg_out.dist.params[:p]
    end
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionMixtureZCCNPX(msg_out::Message{Categorical},
                                       msg_in1::Message{Categorical},
                                       msg_switch::Nothing,
                                       msg_factors::Vararg{Message{PointMass, MatrixVariate}})
    n_factors = length(msg_factors)
    a = zeros(n_factors)
    for k=1:n_factors
        a[k] = msg_out.dist.params[:p]'*msg_factors[k].dist.params[:m]*msg_in1.dist.params[:p]
    end
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSVBTransitionMixtureOutNCCDX(msg_out::Any,
                                          msg_in1::Message{Categorical},
                                          msg_switch::Message{Categorical},
                                          dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})
    z_bar = msg_switch.dist.params[:p]
    n_factors = length(dist_factors)
    d = dims(dist_factors[1])[1] # Dimensionality of output
    a = tiny*ones(d)
    for k=1:n_factors
        a += z_bar[k]*exp.(unsafeLogMean(dist_factors[k]))*msg_in1.dist.params[:p]
    end
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSVBTransitionMixtureIn1CNCDX(msg_out::Message{Categorical},
                                          msg_in1::Any,
                                          msg_switch::Message{Categorical},
                                          dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})
    z_bar = msg_switch.dist.params[:p]
    n_factors = length(dist_factors)
    d = dims(dist_factors[1])[2] # Dimensionality of input
    a = tiny*ones(d)
    for k=1:n_factors
        a += z_bar[k]*exp.(unsafeLogMean(dist_factors[k]))'*msg_out.dist.params[:p]
    end
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSVBTransitionMixtureZCCNDX(msg_out::Message{Categorical},
                                        msg_in1::Message{Categorical},
                                        msg_switch::Any,
                                        dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})
    n_factors = length(dist_factors)
    a = zeros(n_factors)
    for k=1:n_factors
        a[k] = msg_out.dist.params[:p]'*exp.(unsafeLogMean(dist_factors[k]))*msg_in1.dist.params[:p]
    end
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSVBTransitionMixtureA(dist_out_in1_switch::ProbabilityDistribution{Multivariate, Contingency},
                                   dist_factors::Vararg{Union{Nothing, ProbabilityDistribution{MatrixVariate}}})
    
    k = findfirst(dist_factors .== nothing) # Find factor
    
    return Message(MatrixVariate, Dirichlet, a=dist_out_in1_switch.params[:p][k] .+ 1)
end

function ruleMTransitionMixtureCCCDX(msg_out::Message{Categorical, Univariate},
                                     msg_in1::Message{Categorical, Univariate},
                                     msg_switch::Message{Categorical, Univariate},
                                     dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})
    n_factors = length(dist_factors)    
    z_bar = msg_switch.dist.params[:p]
    B = Vector{Matrix}(undef, n_factors)
    for k=1:n_factors
        B[k] = z_bar[k]*Diagonal(msg_out.dist.params[:p])*exp.(unsafeLogMean(dist_factors[k]))*Diagonal(msg_in1.dist.params[:p])
    end

    ProbabilityDistribution(Multivariate, Contingency, p=B./sum(sum(B)))
end

function ruleMTransitionMixtureCCCNX(msg_out::Message{Categorical, Univariate},
                                     msg_in1::Message{Categorical, Univariate},
                                     msg_switch::Message{Categorical, Univariate},
                                     msg_factors::Vararg{Message{PointMass, MatrixVariate}})
    n_factors = length(msg_factors)    
    z_bar = msg_switch.dist.params[:p]
    B = Vector{Matrix}(undef, n_factors)
    for k=1:n_factors
        B[k] = z_bar[k]*Diagonal(msg_out.dist.params[:p])*msg_factors[k].dist.params[:m]*Diagonal(msg_in1.dist.params[:p])
    end

    ProbabilityDistribution(Multivariate, Contingency, p=B./sum(sum(B)))
end
