****************************************
 Probability distributions and messages
****************************************


The ``ProbabilityDistribution`` type
====================================

ForneyLab comes with its own implementation of probability distributions (instead of using `Distributions.jl <https://github.com/JuliaStats/Distributions.jl>`_) for flexibility and efficiency reasons. Every probability distribution type is a subtype of either ``UnivariateProbabilityDistribution`` or ``MultivariateProbabilityDistribution``::

    abstract ProbabilityDistribution
    abstract UnivariateProbabilityDistribution <: ProbabilityDistribution
    abstract MultivariateProbabilityDistribution <: ProbabilityDistribution

The name of a subtype of ``ProbabilityDistribution`` should end in "Distribution", and at least the following methods should be implemented for it: ``==``, ``show``, ``isProper``, and ``vague``. The name of a multivariate distribution type should start with "Mv", i.e. :class:`MvGaussianDistribution`.

.. function:: vague(distribution_type)

    Creates a vague 'almost non-informative' ``ProbabilityDistribution`` of type ``distribution_type``. For the :class:`GaussianDistribution` this means for example a distribution with maximum variance::

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

An :class:`Edge` represents a variable in the probabilistic model. The ``edge.marginal`` field holds the marginal distribution of the corresponding variable. There are some helper functions available to work with marginals.

.. function:: calculateMarginal(edge)

    calculates and returns the marginal distribution from the forward and backward messages present on ``edge``.

.. function:: calculateMarginal!(edge)

    Identical to ``calculateMarginal(edge)``, but the calculated marginal is also written to ``edge.marginal``.
