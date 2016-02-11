type RecognitionDistribution
    distribution::ProbabilityDistribution
    edges::Set{Edge} # Edges on which the distribution is defined
end

function show(io::IO, f::RecognitionDistribution)
    println(io, "RecognitionDistribution: $(f.distribution)")
    println(io, "Defined on edges:")
    for e in f.edges
        println(" $(e)")
    end
end


############################################
# GammaDistribution
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end


############################################
# InverseGammaDistribution
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end


############################################
# BetaDistribution
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end


############################################
# GaussianDistribution
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end


############################################
# MvGaussianDistribution
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::MvGaussianDistribution, backward_dist::MvGaussianDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end


############################################
# WishartDistribution
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::WishartDistribution, backward_dist::WishartDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end


############################################
# Gaussian-students t combination
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end
calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateRecognitionDistribution!(marg, backward_dist, forward_dist)


############################################
# Gaussian-delta combination
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::GaussianDistribution, backward_dist::DeltaDistribution{Float64})
    # Calculation for univariate approximate marginal
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end
calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::DeltaDistribution{Float64}, backward_dist::GaussianDistribution) = calculateRecognitionDistribution!(marg, backward_dist, forward_dist)


############################################
# MvGaussian-MvDelta combination
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::MvGaussianDistribution, backward_dist::MvDeltaDistribution{Float64})
    # Calculation for univariate approximate marginal
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end
calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::MvDeltaDistribution{Float64}, backward_dist::MvGaussianDistribution) = calculateRecognitionDistribution!(marg, backward_dist, forward_dist)


############################
# Joint marginals
############################

function calculateRecognitionDistribution!(q_dist::RecognitionDistribution,
                                 node::GaussianNode,
                                 gaus_msg::Message{GaussianDistribution},
                                 gam_msg::Message{GammaDistribution},
                                 y_dist::GaussianDistribution)
    # (Joint) marginal update function used for SVMP
    # Definitions available in derivations notebook

    marg = q_dist.distribution

    mu_m = gaus_msg.payload
    mu_gam = gam_msg.payload

    ForneyLab.ensureParameters!(mu_m, (:m,))
    ForneyLab.ensureParameters!(y_dist, (:m, :W))

    marg.m = mu_m.m
    marg.beta = huge
    marg.a = mu_gam.a + 0.5
    marg.b = (1.0/(2.0*y_dist.W)) + mu_gam.b

    return marg
end
