****************************************
 Probability distributions and messages
****************************************


The ``ProbabilityDistribution`` type
====================================

ForneyLab comes with its own implementation of probability distributions (instead of using i.e. `Distributions.jl <https://github.com/JuliaStats/Distributions.jl>`_) for flexibility and efficiency reasons. Every probability distribution type is a subtype of either ``UnivariateProbabilityDistribution`` or ``MultivariateProbabilityDistribution``::

    abstract ProbabilityDistribution
    abstract UnivariateProbabilityDistribution <: ProbabilityDistribution
    abstract MultivariateProbabilityDistribution <: ProbabilityDistribution

The name of a subtype of ``ProbabilityDistribution`` should end in "Distribution", and at least the following methods should be implemented for it: ``==``, ``show``, ``isProper``, and ``vague``. The name of a multivariate distribution type should start with "Mv", i.e. :class:`MvGaussianDistribution`.

.. function:: vague(distribution_type)

    Creates a vague 'non-informative' ``ProbabilityDistribution`` of type ``distribution_type``. For the :class:`GaussianDistribution` this means for example a distribution with maximum variance::

        non_informative_gaussian = vague(GaussianDistribution) # identical to GaussianDistribution(m=0.0, V=huge)

.. function:: isProper(p::ProbabilityDistribution)

    Returns true/false according to whether or not ``p`` is a proper probability distribution. A distribution is proper if and only if (i) the pdf/pmf is upper bounded by 1, and (ii) the integral/sum of the pdf/pmf over the entire domain equals 1.


Moreover, the following optional methods might be implemented:

.. function:: mean(p::ProbabilityDistribution)

    Returns the mean of ``p``.


.. function:: var(p::ProbabilityDistribution)

    Returns the variance of ``p``.


.. function:: sample(p::ProbabilityDistribution)

    Returns a random sample drawn from ``p``.


Built-in distributions
======================

Univariate distributions
------------------------

Built-in univariate distributions: :class:`BetaDistribution`, :class:`DeltaDistribution`, :class:`GammaDistribution`, :class:`GaussianDistribution`, :class:`InverseGammaDistribution`, :class:`StudentsTDistribution`.


.. type:: BetaDistribution

    :description:   Beta distribution (univariate)
    :parameters:    ``a > 0`` ("shape", real scalar), ``b > 0`` ("rate", real scalar)
    :construction:  ``BetaDistribution(a=1.0, b=1.0)``
    :reference:     Bishop, 2006; Pattern recognition and machine learning; appendix B

.. type:: DeltaDistribution

    :description:   Kronecker delta (``pdf(x) = δ(x-m)``)
    :parameters:    ``m`` (Any)
    :construction:  ``DeltaDistribution(m)``

    The ``DeltaDistribution`` is used to fix variables to a value, for example to capture observed data.


.. type:: GammaDistribution

    :description:   Gamma distribution (univariate)
    :parameters:    ``a > 0`` ("shape", real scalar), ``b > 0`` ("rate", real scalar)
    :construction:  ``GammaDistribution(a=1.0, b=1.0)``
    :reference:     Bishop, 2006; Pattern recognition and machine learning; appendix B

.. type:: LogNormalDistribution

    :description:   Log-normal distribution (univariate)
    :parameters:    ``m`` ("location", real scalar), ``s > 0`` ("squared scale" (s = σ²), real scalar)
    :construction:  ``LogNormalDistribution(m=0.0, s=1.0)``

.. type:: GaussianDistribution

    :description:   Gaussian distribution (multivariate)
    :parameters:    ``m`` ("mean", real scalar), ``V`` ("variance", real scalar), ``W`` ("precision", real scalar), ``xi`` ("weighted mean", real scalar)
    :construction:  ``GaussianDistribution(m=0.0, V=1.0)`` or ``GaussianDistribution(xi=0.0, W=1.0)`` or any other valid parameter combination.
    :reference:     Bishop, 2006; Pattern recognition and machine learning; appendix B

    The Gaussian distribution can be parametrized in multiple ways. Depending on the application, a specific parametrization might be attractive from a computational point of view. The following combinations are valid: ``(m,V)``, ``(m,W)``, ``(xi,V)``, ``(xi,W)``.
    As long as ``V`` and ``W`` are non-zero, the parametrizations can be converted using:

    .. math::
        \begin{aligned}
        W &= V^{-1} \\
        ξ &= W⋅m
        \end{aligned}

    The following functions are available to facilitate parameter conversions:

    .. function:: ensureParameters!(dist::GaussianDistribution, params::Tuple{Symbol})

        Make sure that the specified parameters of ``dist`` are set and valid. Calculate them from the other (valid) parameters if required.
        Example: ``ensureParameters!(dist, (:m,:V))`` to make sure that ``dist.m`` and ``dist.V`` are set.

    .. function:: isWellDefined(dist::GaussianDistribution)

        Returns ``true`` if and only if ``dist`` has a valid parametrization.

    .. function:: isConsistent(dist::GaussianDistribution)

        Returns ``true`` if and only if ``dist`` has a consistent parametrization. If ``dist`` is 'overdetermined', this function validates that the multiple parametrizations are in agreement.

    The parameters of ``GaussianDistribution`` that are *not* used or that are not valid should be invalidated by setting them to ``NaN``. Validity of a parameter should be checked using ``isNaN(...)``.


.. type:: InverseGammaDistribution

    :description:   Inverse-gamma distribution (univariate)
    :parameters:    ``a > 0`` ("shape", real scalar), ``b > 0`` ("scale", real scalar)
    :construction:  ``InverseGammaDistribution(a=1.0, b=1.0)``
    :reference:     Korl, 2005; A factor graph approach to signal modelling, system identification and filtering; appendix A


.. type:: StudentsTDistribution

    :description:   Student's t-distribution (multivariate)
    :parameters:    ``m`` ("mean", real vector), ``lambda`` ("inverse scale", positive definite real matrix), ``nu`` ("degrees of freedom", real scalar)
    :construction:  ``StudentsTDistribution(m, lambda, nu)``
    :reference:     Bishop, 2006; Pattern recognition and machine learning; appendix B


Multivariate distributions
--------------------------

Built-in multivariate distributions: :class:`MvDeltaDistribution`, :class:`MvGaussianDistribution`, :class:`MvNormalGammaDistribution`.

.. type:: MvDeltaDistribution

    :description:   Multivariate Kronecker delta (``pdf(x) = δ(x-m)``)
    :parameters:    ``m`` (``Vector{Any}``)
    :construction:  ``MvDeltaDistribution(m)``

    Same as :class:`DeltaDistribution`, just for the multivariate case.


.. type:: MvGaussianDistribution

    :description:   Gaussian distribution (multivariate)
    :parameters:    ``m`` ("mean", real vector), ``V`` ("variance", real matrix), ``W`` ("precision", real matrix), ``xi`` ("weighted mean", real vector)
    :construction:  ``MvGaussianDistribution(m=zeros(3), V=eye(3))`` or ``MvGaussianDistribution(xi=zeros(3), W=2.0*eye(3))`` or any other valid parameter combination.
    :reference:     Bishop, 2006; Pattern recognition and machine learning; appendix B

    The parametrization options for this distribution are equal to those of :class:`GaussianDistribution`.

    The following functions are available to facilitate parameter conversions:

    .. function:: ensureParameters!(dist::MvGaussianDistribution, params::Tuple{Symbol})

        Make sure that the specified parameters of ``dist`` are set and valid. Calculate them from the other (valid) parameters if required.
        Example: ``ensureParameters!(dist, (:m,:V))`` to make sure that ``dist.m`` and ``dist.V`` are set.

    .. function:: isWellDefined(dist::MvGaussianDistribution)

        Returns ``true`` if and only if ``dist`` has a valid parametrization.

    .. function:: isConsistent(dist::MvGaussianDistribution)

        Returns ``true`` if and only if ``dist`` has a consistent parametrization. If ``dist`` is 'overdetermined', this function validates that the multiple parametrizations are in agreement.

    The parameters of ``MvGaussianDistribution`` that are *not* used or that are not valid should be invalidated using :func:`invalidate!()`. Validity of a parameter can be checked using :func:`isValid()`.


.. type:: MvNormalGammaDistribution

    :description:   Normal-gamma distribution (bivariate)
    :parameters:    ``m`` ("location", real scalar), ``beta > 0`` ("precision", real scalar), ``a`` ("shape", real scalar), ``b`` ("rate", real scalar)
    :construction:  ``MvNormalGammaDistribution(m=0.0, beta=1.0, a=1.0, b=1.0)``
    :reference:     Bishop, 2006; Pattern recognition and machine learning; appendix B



Messages
========

.. type:: Message

    ::

        type Message{T<:ProbabilityDistribution}
            payload::T
        end

    Messages are passed over edges, and carry a :class:`ProbabilityDistribution` in the ``payload`` field. A ``Message`` is usually stored on an :class:`Interface`.


Marginals
=========

.. seealso::
    **Demo:** `Marginals <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/02_marginals.ipynb>`_

Since an :class:`Edge` represents a variable in the probabilistic model, the ``edge.marginal`` field holds the marginal distribution of the corresponding variable. There are some helper functions available to work with marginals.

.. function:: calculateMarginal(edge)

    If the forward and backward messages on ``edge`` are calculated according to the sum-product rule, the marginal distribution of the variable represented by ``edge`` can be calculated from these messages. ``calculateMarginal(edge)`` calculates and returns the marginal distribution from the forward and backward messages.

.. function:: calculateMarginal!(edge)

    Identical to ``calculateMarginal(edge)``, but the calculated marginal is also written to ``edge.marginal``.


.. function:: getMarginalType(distributions...)

    Returns the type of the marginal distribution given the types of its factors (i.e. carried by forward/backward messages).
