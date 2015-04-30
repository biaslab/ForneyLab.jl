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
    return ForneyLab.gammaMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# InverseGammaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    return ForneyLab.inverseGammaMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# BetaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    return ForneyLab.betaMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# GaussianDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    return ForneyLab.gaussianMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# Gaussian-students t combination
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    return ForneyLab.gaussianStudentsMarginalRule!(marg.distribution, forward_dist, backward_dist)
end
calculateQDistribution!(marg::QDistribution, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(factor, graph, backward_dist, forward_dist)

# Gaussian-Delta combination
function calculateQDistribution!(marg::QDistribution, forward_dist::DeltaDistribution, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    return marg.distribution = convert(GaussianDistribution, forward_dist)
end
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::DeltaDistribution)
    # Calculation for univariate approximate marginal
    return marg.distribution = convert(GaussianDistribution, backward_dist)
end


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