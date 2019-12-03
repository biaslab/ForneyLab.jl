export ruleVBLogitOut, ruleVBLogitIn1, ruleVBLogitXi

function ruleVBLogitOut(marg_out::Any, 
                        marg_in1::ProbabilityDistribution{Univariate}, 
                        marg_xi::ProbabilityDistribution{Univariate})

    return Message(Univariate, Bernoulli, p=logisticSigmoid(unsafeMean(marg_in1)))
end

function ruleVBLogitIn1(marg_out::ProbabilityDistribution{Univariate}, 
                        marg_in1::Any, 
                        marg_xi::ProbabilityDistribution{Univariate})
    
    xi_hat = unsafeMode(marg_xi)
    
    return Message(Univariate, GaussianWeightedMeanPrecision, xi=unsafeMean(marg_out) - 0.5, w=2*logisticLambda(xi_hat))
end

function ruleVBLogitXi(marg_out::ProbabilityDistribution{Univariate}, 
                       marg_in1::ProbabilityDistribution{Univariate}, 
                       marg_xi::Any)
    
    return Message(Univariate, Function, mode=sqrt(unsafeMean(marg_in1)^2 + unsafeCov(marg_in1)))
end