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

function calculateRecognitionDistribution!(recog::RecognitionDistribution, forward_dist::ProbabilityDistribution, backward_dist::ProbabilityDistribution)
    # Recognition distribution is the product of the forward and backward distributions
    return prod!(forward_dist, backward_dist, recog.distribution)
end

############################################
# Bernoulli-Bernoulli combination
############################################

function calculateRecognitionDistribution!(marg::RecognitionDistribution, forward_dist::BernoulliDistribution, backward_dist::BernoulliDistribution)
    # Calculation for univariate approximate marginal
    return ForneyLab.equalityRule!(marg.distribution, forward_dist, backward_dist)
end

############################
# Joint marginals
############################

function calculateRecognitionDistribution!(q_dist::RecognitionDistribution,
                                 node::GaussianNode,
                                 gaus_msg::Message{Gaussian},
                                 gam_msg::Message{Gamma},
                                 y_dist::Gaussian)
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
