export ruleVBLogitOut, ruleVBLogitIn1, ruleVBLogitXi

function ruleVBLogitOut(marg_out::Any, 
                        marg_in1::Distribution{Univariate}, 
                        marg_xi::Distribution{Univariate})

    return Message(Univariate, Bernoulli, p=logisticSigmoid(unsafeMean(marg_in1)))
end

function ruleVBLogitIn1(marg_out::Distribution{Univariate}, 
                        marg_in1::Any, 
                        marg_xi::Distribution{Univariate})
    
    xi_hat = unsafeMode(marg_xi)
    
    return Message(Univariate, Gaussian{Canonical}, xi=unsafeMean(marg_out) - 0.5, w=2*logisticLambda(xi_hat))
end

function ruleVBLogitXi(marg_out::Distribution{Univariate}, 
                       marg_in1::Distribution{Univariate}, 
                       marg_xi::Any)
    
    return Message(Univariate, Function, mode=sqrt(unsafeMean(marg_in1)^2 + unsafeCov(marg_in1)))
end