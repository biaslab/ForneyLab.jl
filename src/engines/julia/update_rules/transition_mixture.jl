export
ruleSPTransitionMixtureOutNCPX,
ruleSPTransitionMixtureIn1CNPX,
ruleVBTransitionMixtureZ, 
ruleVBTransitionMixtureOut,
ruleVBTransitionMixtureIn1,
ruleVBTransitionMixtureA

function ruleSPTransitionMixtureOutNCPX(msg_out::Nothing,
                                        msg_in1::Message{Categorical},
                                        msg_switch::Message{PointMass},
                                        msg_factors::Vararg{Message{PointMass, MatrixVariate}})
    
    i = findfirst(Bool.(unsafeMeanVector(msg_switch.dist)))
    a = clamp.(msg_factors[i].dist.params[:m]*msg_in1.dist.params[:p], tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleSPTransitionMixtureIn1CNPX(msg_out::Message{Categorical},
                                        msg_in1::Nothing,
                                        msg_switch::Message{PointMass},
                                        msg_factors::Vararg{Message{PointMass, MatrixVariate}})
    
    i = findfirst(Bool.(unsafeMeanVector(msg_switch.dist)))
    a = clamp.(msg_factors[i].dist.params[:m]'*msg_out.dist.params[:p], tiny, Inf)
    
    Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionMixtureZ(dist_out::ProbabilityDistribution,
                                  dist_in1::ProbabilityDistribution,
                                  dist_switch::Any,
                                  dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})

    n_factors = length(dist_factors)
    U = Vector{Float64}(undef, n_factors)
    for k = 1:n_factors
        U[k] = averageEnergy(Transition, dist_out, dist_in1, dist_factors[k])
    end

    return Message(Univariate, Categorical, p=softmax(-U))
end

function ruleVBTransitionMixtureOut(dist_out::Any,
                                    dist_in1::ProbabilityDistribution,
                                    dist_switch::ProbabilityDistribution,
                                    dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})

    z_bar = unsafeMeanVector(dist_switch)
    d = dims(dist_factors[1])[1] # Dimensionality of output
    log_a = zeros(d)
    for k = 1:length(z_bar)
        log_a += z_bar[k]*unsafeLogMean(dist_factors[k])*unsafeMeanVector(dist_in1)
    end

    a = clamp.(exp.(log_a), tiny, Inf)
    
    return Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionMixtureIn1(dist_out::ProbabilityDistribution,
                                    dist_in1::Any,
                                    dist_switch::ProbabilityDistribution,
                                    dist_factors::Vararg{ProbabilityDistribution{MatrixVariate}})

    z_bar = unsafeMeanVector(dist_switch)
    d = dims(dist_factors[1])[2] # Dimensionality of input
    log_a = zeros(d)
    for k = 1:length(z_bar)
        log_a += z_bar[k]*unsafeLogMean(dist_factors[k])'*unsafeMeanVector(dist_out)
    end

    a = clamp.(exp.(log_a), tiny, Inf)
    
    return Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBTransitionMixtureA(dist_out::ProbabilityDistribution,
                                  dist_in1::ProbabilityDistribution,
                                  dist_switch::ProbabilityDistribution,
                                  dist_factors::Vararg{Union{Nothing, ProbabilityDistribution{MatrixVariate}}})
    
    k = findfirst(dist_factors .== nothing) # Find factor
    z_bar = unsafeMeanVector(dist_switch)
    
    return Message(MatrixVariate, Dirichlet, a=z_bar[k]*unsafeMeanVector(dist_out)*unsafeMeanVector(dist_in1)' .+ 1)
end


# TODO: Structured updates