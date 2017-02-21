export PointMass

"""Encodes a probability distribution as a FactorNode of type `family` with fixed interfaces"""
immutable ProbabilityDistribution{family<:FactorNode}
    parameters::Tuple
end

mean(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
var(dist::ProbabilityDistribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::ProbabilityDistribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")

"""
PointMass is an abstract node type used to describe point mass distributions.
It never occurs in a FactorGraph, but it is used as a probability distribution type.
"""
abstract PointMass <: DeltaFactor

mean(dist::ProbabilityDistribution{PointMass}) = dist.parameters[1]

isProper(::ProbabilityDistribution{PointMass}) = true