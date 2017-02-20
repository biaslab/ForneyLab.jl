export PointMass

"""
Description:

    Defines a point mass:

    f(x,m) = δ(x - m)

Interfaces:

    1. i[1], 2. i[2]

Construction:

    No constructor available.
"""
abstract PointMass <: DeltaFactor

mean(dist::ProbabilityDistribution{PointMass}) = dist.parameters[1]

isProper(::ProbabilityDistribution{PointMass}) = true