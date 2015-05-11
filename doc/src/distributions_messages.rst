****************************************
 Probability distributions and messages
****************************************


The ``ProbabilityDistribution`` type
====================================

ForneyLab comes with its own implementation of probability distributions (instead of using i.e. `Distributions.jl <https://github.com/JuliaStats/Distributions.jl>`_) for flexibility and efficiency reasons. Every probability distribution type is a subtype of::

    abstract ProbabilityDistribution

The name of a subtype of ``ProbabilityDistribution`` should end in "Distribution", and the following methods should be implemented for it: ``==``, ``show``, and ``vague``. 

.. function:: vague(distribution_type)

    Creates a vague 'non-informative' ``ProbabilityDistribution`` of type ``distribution_type``. For the :class:`GaussianDistribution` this means for example a distribution with maximum variance::

        non_informative_gaussian = vague(GaussianDistribution) # identical to GaussianDistribution(m=0.0, V=huge())

Moreover, the following optional methods might be implemented:

.. function:: mean(p::ProbabilityDistribution)

    Returns the first moment of ``p``.


.. function:: var(p::ProbabilityDistribution)

    Returns the second moment (variance) of ``p``.


.. function:: sample(p::ProbabilityDistribution)

    Returns a random sample drawn from ``p``.


Built-in distributions
======================

The following built-in probability distributions are available: :class:`BetaDistribution`, :class:`DeltaDistribution`, :class:`GammaDistribution`, :class:`GaussianDistribution`, :class:`InverseGammaDistribution`, :class:`NormalGammaDistribution`, :class:`StudentsTDistribution`.

.. type:: BetaDistribution

    :description:   Beta distribution (univariate)
    :parameters:    ``a > 0`` ("shape", Real scalar), ``b > 0`` ("rate", Real scalar)
    :construction:  ``BetaDistribution(a=1.0, b=1.0)``


.. type:: DeltaDistribution

    :description:   Dirac/Kronecker delta (``pdf(x) = δ(x-m)``)
    :parameters:    ``m`` (Any)
    :construction:  ``DeltaDistribution(m)``

    The ``DeltaDistribution`` is used to fix variables to a value, for example to capture observed data.


.. type:: GammaDistribution

    :description:   Gamma distribution (univariate)
    :parameters:    ``a > 0`` ("shape", Real scalar), ``b > 0`` ("rate", Real scalar)
    :construction:  ``GammaDistribution(a=1.0, b=1.0)``


.. type:: GaussianDistribution

    :description:   Gaussian distribution (multivariate)
    :parameters:    ``m`` ("mean", Real vector), ``V`` ("variance", Real positive definite matrix), ``W`` ("precision", Real positive definite matrix), ``xi`` ("weighted mean", Real vector)
    :construction:  ``GaussianDistribution(m=0.0, V=1.0)`` or ``GaussianDistribution(xi=0.0, W=1.0)`` or any other valid parameter combination.

    The Gaussian distribution can be parametrized in multiple ways. Depending on the application, a specific parametrization might be attractive from a computational point of view. The following combinations are valid: ``(m,V)``, ``(m,W)``, ``(xi,V)``, ``(xi,W)``.
    As long as ``V`` and ``W`` are non-singular, the parametrizations can be converted using:

    .. math:: 
        \begin{aligned}
        W &= V^{-1} \\
        ξ &= W⋅m
        \end{aligned}

    The following functions are available to facilitate parameter conversions:

    .. function:: ensureMVParametrization!(dist::GaussianDistribution)

        Make sure ``dist.m`` and ``dist.V`` are defined and valid. Calculate from other parameters if required.

    .. function:: ensureMWParametrization!(dist::GaussianDistribution)

        Make sure ``dist.m`` and ``dist.W`` are defined and valid. Calculate from other parameters if required.

    .. function:: ensureXiVParametrization!(dist::GaussianDistribution)

        Make sure ``dist.xi`` and ``dist.V`` are defined and valid. Calculate from other parameters if required.

    .. function:: ensureXiWParametrization!(dist::GaussianDistribution)

        Make sure ``dist.xi`` and ``dist.W`` are defined and valid. Calculate from other parameters if required.

    .. function:: isWellDefined(dist::GaussianDistribution)

        Returns ``true`` if and only if ``dist`` has a valid parametrization.

    .. function:: isConsistent(dist::GaussianDistribution)

        Returns ``true`` if and only if ``dist`` has a consistent parametrization. If ``dist`` is 'overdetermined', this function validates that the multiple parametrizations are in agreement.

    The parameters of ``GaussianDistribution`` that are *not* used or that are not valid should be invalidated using :func:`invalidate!()`. Validity of a parameter can be checked using :func:`isValid()`.


.. type:: InverseGammaDistribution

    :description:   Inverse-gamma distribution (univariate)
    :parameters:    ``a > 0`` ("shape", Real scalar), ``b > 0`` ("scale", Real scalar)
    :construction:  ``InverseGammaDistribution(a=1.0, b=1.0)``


.. type:: NormalGammaDistribution

    :description:   Normal-gamma distribution (bivariate)
    :parameters:    ``m`` ("location", Real scalar), ``beta > 0`` ("precision", Real scalar), ``a`` ("shape", Real scalar), ``b`` ("rate", Real scalar)
    :construction:  ``NormalGammaDistribution(m=0.0, beta=1.0, a=1.0, b=1.0)``


.. type:: StudentsTDistribution

    :description:   Student's t-distribution (multivariate)
    :parameters:    ``m`` ("mean", Real vector), ``W`` ("precision", positive definite Real matrix), ``nu`` ("degrees of freedom", Real scalar)
    :construction:  ``StudentsTDistribution(m, W, nu)``


Messages
========

.. type:: Message



Marginals
=========

.. function:: calculateMarginal(edge)


.. function:: calculateMarginal!(edge)


.. function:: getMarginalType(distributions...)

