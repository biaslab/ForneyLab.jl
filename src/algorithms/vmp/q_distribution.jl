type QDistribution
    distribution::ProbabilityDistribution
    edges::Set{Edge} # Edges on which the distribution is defined
end

function show(io::IO, f::QDistribution)
    println(io, "QDistribution: $(f.distribution)")
    println(io, "Defined on edges:")
    for e in f.edges
        println(" $(e)")
    end
end


############################
# Univariate approximate marginal calculations
############################

# GammaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    return ForneyLab.equalityGammaRule!(marg.distribution, forward_dist, backward_dist)
end

# InverseGammaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    return ForneyLab.equalityInverseGammaRule!(marg.distribution, forward_dist, backward_dist)
end

# BetaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    return ForneyLab.equalityBetaRule!(marg.distribution, forward_dist, backward_dist)
end

# GaussianDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    return ForneyLab.equalityGaussianRule!(marg.distribution, forward_dist, backward_dist)
end

# Gaussian-students t combination
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    return ForneyLab.equalityGaussianStudentsTRule!(marg.distribution, forward_dist, backward_dist)
end
calculateQDistribution!(marg::QDistribution, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateQDistribution!(marg, backward_dist, forward_dist)

# Gaussian-delta combination
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::DeltaDistribution)
    # Calculation for univariate approximate marginal
    return marg.distribution = ForneyLab.equalityGaussianDeltaRule!(marg.distribution, forward_dist, backward_dist)
end
calculateQDistribution!(marg::QDistribution, forward_dist::DeltaDistribution, backward_dist::GaussianDistribution) = calculateQDistribution!(marg, backward_dist, forward_dist)

# TODO: check for need
# # Gamma-delta combination
# function calculateQDistribution!(marg::QDistribution, forward_dist::GammaDistribution, backward_dist::DeltaDistribution)
#     # Calculation for univariate approximate marginal
#     return marg.distribution = ForneyLab.equalityGammaDeltaRule!(marg.distribution, forward_dist, backward_dist)
# end
# calculateQDistribution!(marg::QDistribution, forward_dist::DeltaDistribution, backward_dist::GammaDistribution) = calculateQDistribution!(marg, backward_dist, forward_dist)

# TODO: check for need
# # Inverse gamma-delta combination
# function calculateQDistribution!(marg::QDistribution, forward_dist::InverseGammaDistribution, backward_dist::DeltaDistribution)
#     # Calculation for univariate approximate marginal
#     return marg.distribution = ForneyLab.equalityInverseGammaDeltaRule!(marg.distribution, forward_dist, backward_dist)
# end
# calculateQDistribution!(marg::QDistribution, forward_dist::DeltaDistribution, backward_dist::InverseGammaDistribution) = calculateQDistribution!(marg, backward_dist, forward_dist)

############################
# Joint approximate marginal calculations
############################

function calculateQDistribution!(q_dist::QDistribution,
                                 node::GaussianNode,
                                 gaus_msg::Message{GaussianDistribution},
                                 gam_msg::Message{GammaDistribution},
                                 y_dist::GaussianDistribution)
    # (Joint) marginal update function used for SVMP
    # Definitions available in derivations notebook

    marg = q_dist.distribution

    mu_m = gaus_msg.payload
    mu_gam = gam_msg.payload
    
    ForneyLab.ensureMDefined!(mu_m)
    ForneyLab.ensureMWParametrization!(y_dist)

    (length(mu_m.m) == 1 && length(y_dist.m) == 1)|| error("Update rule for NormalGammaDistribution marginal only supports univariate distributions.")

    marg.m = mu_m.m[1]
    marg.beta = huge()
    marg.a = mu_gam.a + 0.5
    marg.b = (1.0/(2.0*y_dist.W[1, 1])) + mu_gam.b

    return marg
end