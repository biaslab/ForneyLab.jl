export ruleVBSoftmaxOut, ruleVBSoftmaxIn1, ruleVBSoftmaxXi, ruleVBSoftmaxA

function ruleVBSoftmaxOut(marg_out::Any, 
                          marg_in1::ProbabilityDistribution{Multivariate}, 
                          marg_xi::ProbabilityDistribution{Multivariate, PointMass},
                          marg_a::ProbabilityDistribution{Univariate, PointMass})
    
    b_bar = unsafeBoundMean(marg_in1, marg_xi, marg_a)
    a = exp.(unsafeMean(marg_in1) .- b_bar)

    return Message(Univariate, Categorical, p=a./sum(a))
end

function ruleVBSoftmaxIn1(marg_out::ProbabilityDistribution{Univariate}, 
                          marg_in1::Any, 
                          marg_xi::ProbabilityDistribution{Multivariate, PointMass},
                          marg_a::ProbabilityDistribution{Univariate, PointMass})
    
    xi_hat = marg_xi.params[:m]
    a_hat = marg_a.params[:m]
    lambda_xi_hat = logisticLambda.(xi_hat)

    gam = 2*unsafeMean(marg_out)*lambda_xi_hat'
    mu_plus = a_hat .+ (1.0./(4*lambda_xi_hat))
    mu_min = a_hat .- (1.0./(4*lambda_xi_hat))

    xi = mu_plus.*diag(gam) + (gam - Diagonal(diag(gam)))'*mu_min
    W = Diagonal(vec(sum(gam, dims=1)))

    return Message(Multivariate, GaussianWeightedMeanPrecision, xi=xi, w=W)
end

function ruleVBSoftmaxXi(marg_out::ProbabilityDistribution{Univariate}, 
                         marg_in1::ProbabilityDistribution{Multivariate}, 
                         marg_xi::Any,
                         marg_a::ProbabilityDistribution{Univariate, PointMass})
    
    return Message(Multivariate, PointMass, m=sqrt.(unsafeVar(marg_in1) + (unsafeMean(marg_in1) .- marg_a.params[:m]).^2))
end

function ruleVBSoftmaxA(marg_out::ProbabilityDistribution{Univariate}, 
                        marg_in1::ProbabilityDistribution{Multivariate}, 
                        marg_xi::ProbabilityDistribution{Multivariate, PointMass},
                        marg_a::Any)
    
    xi_hat = marg_xi.params[:m]
    lambda_xi_hat = logisticLambda.(xi_hat)
    d = length(xi_hat)

    return Message(Univariate, PointMass, m=((d-2)/4 + unsafeMean(marg_in1)'*lambda_xi_hat)/sum(lambda_xi_hat))
end