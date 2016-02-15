**************
 Introduction
**************

Several pages in the documentation are accompanied with demos that showcase ForneyLab functionality in practice. Relevant demos are referenced at the beginning of a section. For example, an introductory demo is available here:

.. seealso::
    **Demo:** `Basics <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/01_basics.ipynb>`_

A factor graph is a graphical model that encodes a factorization of the joint probability distribution over all variables involved. In a Forney-style factor graph, every variable is represented by an edge. The nodes connecting the edges represent the factors.

We will introduce the basic building blocks using a simple example. Consider the following probabilistical model (random walk):

.. math::
    \begin{aligned}
    X_2 &= X_1 + \epsilon_1, \\
    X_3 &= X_2 + \epsilon_2, \\
    \epsilon_1 &\sim \mathcal{N}(0,1),\\
    \epsilon_2 &\sim \mathcal{N}(2,1).
    \end{aligned}

This model defines a factorization of the joint probability distribution:

.. math::
    p(X_1,X_2,X_3,\epsilon_1,\epsilon_2) = p(X_1,X_2,\epsilon_1) \, p(X_2,X_3,\epsilon_2) \, p(X1) \, p(X2) \, p(X3) \, p(\epsilon_1) \, p(\epsilon_2).

The corresponding FFG is::

          | ε1    | ε2
          |       |
      X1  v   X2  v   X3
    ---->[+]---->[+]---->

The graph is directed to indicate direction of the "+" operation. Note that only X2 is a full edge, the others are half-edges since they are connected to just one node. The general flow when working with FFGs is:

1. Building a FFG to model a specific system.
2. Fixing some variables in the model and deriving an algorithm for inferring the values of the other variables given the fixed ones.
3. Running the derived algorithm.

In the case of our running example, we might be interested in inferring the value of X3 given the noise distributions and (an estimate of) of X1. Inference has different meanings in different contexts:

- Parameter estimation / model fitting / system identification: inferring parameter values from data.
- State estimation: inferring (hidden) state variables of a system from data.
- Prediction: inferring future data (from past data and/or parameter values).
- Smoothing: inferring data points from other data and/or parameters.

Inference in a FFG can mean any of the above, depending on which variables (edges) are fixed and which ones are inferred. Inference in a FFG can often be implemented efficiently by a message passing algorithm. ForneyLab provides a framework to:

1. Build (hierarchical) FFGs;
2. Automatically generate message passing algorithms to solve inference problems;
3. Efficiently run message passing algorithms on FFGs.