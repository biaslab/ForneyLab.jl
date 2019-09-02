ForneyLab.jl
============

*A Julia toolbox for automatic generation of (Bayesian) inference algorithms.*

Given a probabilistic model, ForneyLab generates efficient Julia code for message-passing based inference. It uses the model structure to generate an algorithm that consists of a sequence of local computations on a Forney-style factor graph (FFG) representation of the model.

## Package Features

- User friendly syntax for specification of probabilistic models.
- Automatic generation of message passing algorithms including
    - [Belief propagation](https://en.wikipedia.org/wiki/Belief_propagation)
    - [Variational message passing](https://en.wikipedia.org/wiki/Variational_message_passing)
    - [Expectation maximization](https://en.wikipedia.org/wiki/Expectation-maximization_algorithm)
    - [Expectation propagation](https://en.wikipedia.org/wiki/Expectation_propagation)
- Support for hybrid models combining discrete and continuous latent variables.
- Evaluation of free energy as a model performance measure.
- Combination of distinct inference algorithms under a unified paradigm.
- Features composite nodes that allow for flexible hierarchical design in terms of model structure and algorithms.

## Resources

- For an in depth overview of ForneyLab, see [A Factor Graph Approach to Automated Design of Bayesian Signal Processing Algorithms](https://arxiv.org/abs/1811.03407) by Cox et. al. (2018).
- For an introduction to message passing and FFGs, see [The Factor Graph Approach to Model-Based Signal Processing](https://ieeexplore.ieee.org/document/4282128/) by Loeliger et al. (2007).
- The [ForneyLab project page](http://forneylab.org) provides more background on ForneyLab as well as pointers to related literature and talks.

## How to get started?
Head to the [Getting started](@ref) section to get up and running with ForneyLab in no time.
