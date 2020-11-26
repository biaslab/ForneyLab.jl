export
ruleVBTransitionMixtureZ, 
ruleVBTransitionMixtureOut


function ruleVBTransitionMixtureZ(dist_out::ProbabilityDistribution,
                                  dist_switch::Any,
                                  dist_factors::Vararg{ProbabilityDistribution{MatrixVariate, PointMass}})

    n_factors = length(dist_factors)
    U = Vector{Float64}(undef, n_factors)
    for k = 1:n_factors
        A_k = clamp.(dist_factors.params[:m], tiny, 1-tiny) # Soften given transition functions
        A_k = A_k./sum(A_k,dims=1) # Re-normalize columns
        U[k] = averageEnergy(Dirichlet, dist_out, A_k)
    end

    return Message(Univariate, Categorical, p=softmax(-U))
end

function ruleVBTransitionMixtureOut(dist_out::Any,
                                    dist_switch::ProbabilityDistribution,
                                    dist_factors::Vararg{ProbabilityDistribution{MatrixVariate, PointMass}})

    z_bar = unsafeMeanVector(dist_switch)
    A_bar = [dist_factor.params[:m] for dist_factor in dist_factors]

    return Message(MatrixVariate, Dirichlet, a=sum(z_bar.*A_bar))
end