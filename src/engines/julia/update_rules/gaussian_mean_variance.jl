export
ruleSPGaussianMeanVarianceOutPP,
ruleSPGaussianMeanVarianceMPP,
ruleSPGaussianMeanVarianceOutGP, 
ruleSPGaussianMeanVarianceMPG, 
ruleVBGaussianMeanVarianceM,
ruleVBGaussianMeanVarianceOut

ruleSPGaussianMeanVarianceOutPP(msg_out::Void,
                                msg_mean::Message{Univariate{PointMass}},
                                msg_var::Message{Univariate{PointMass}}) =
    Message(Univariate(Gaussian, m=msg_mean.dist.params[:m], v=msg_var.dist.params[:m]))

ruleSPGaussianMeanVarianceMPP(msg_out::Message{Univariate{PointMass}}, msg_mean::Void, msg_var::Message{Univariate{PointMass}}) = ruleSPGaussianMeanVarianceOutPP(msg_mean, msg_out, msg_var)

function ruleSPGaussianMeanVarianceOutGP(   msg_out::Void,
                                            msg_mean::Message{Univariate{Gaussian}},
                                            msg_var::Message{Univariate{PointMass}})

    ensureParameters!(msg_mean.dist, (:m, :v))

    Message(Univariate(Gaussian, m=msg_mean.dist.params[:m], v=msg_mean.dist.params[:v] + msg_var.dist.params[:m]))
end

ruleSPGaussianMeanVarianceMPG(msg_out::Message{Univariate{Gaussian}}, msg_mean::Void, msg_var::Message{Univariate{PointMass}}) = ruleSPGaussianMeanVarianceOutGP(msg_mean, msg_out, msg_var)

ruleVBGaussianMeanVarianceM(dist_out::Univariate,
                            dist_mean::Any,
                            dist_var::Univariate) =
    Message(Univariate(Gaussian, m=unsafeMean(dist_out), v=unsafeMean(dist_var)))

ruleVBGaussianMeanVarianceOut(  dist_out::Any,
                                dist_mean::Univariate,
                                dist_var::Univariate) =
    Message(Univariate(Gaussian, m=unsafeMean(dist_mean), v=unsafeMean(dist_var)))
