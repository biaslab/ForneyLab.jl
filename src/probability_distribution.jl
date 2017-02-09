"""Encodes a probability distribution as a FactorNode of type `family` with fixed interfaces"""
immutable ProbabilityDistribution{family<:FactorNode}
    parameters::Tuple
end

mean(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
var(dist::ProbabilityDistribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::ProbabilityDistribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")