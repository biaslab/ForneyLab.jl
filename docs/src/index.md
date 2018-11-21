ForneyLab.jl
============

[![Build Status](https://travis-ci.org/biaslab/ForneyLab.jl.svg?branch=master)](https://travis-ci.org/biaslab/ForneyLab.jl)

ForneyLab.jl is a Julia package for automatic generation of (Bayesian) inference algorithms. Given a probabilistic model, ForneyLab generates efficient Julia code for message-passing based inference. It uses the model structure to generate an algorithm that consists of a sequence of local computations on a Forney-style factor graph (FFG) representation of the model. For an excellent introduction to message passing and FFGs, see [The Factor Graph Approach to Model-Based Signal Processing](https://ieeexplore.ieee.org/document/4282128/) by Loeliger et al. (2007).

We designed ForneyLab with a focus on flexible and modular modeling of time-series data. ForneyLab enables a user to:

- Conveniently specify a probabilistic model;
- Automatically generate an efficient inference algorithm;
- Compile the inference algorithm to executable Julia code.

The current version supports [belief propagation](https://en.wikipedia.org/wiki/Belief_propagation) (sum-product message passing), [variational message passing](https://en.wikipedia.org/wiki/Variational_message_passing) and [expectation propagation](https://en.wikipedia.org/wiki/Expectation_propagation).

The [ForneyLab project page](http://forneylab.org) provides more background on ForneyLab as well as pointers to related literature and talks. For a practical introduction, have a look at [the demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo).


## Installation

Install ForneyLab through the Julia package manager:
```jl
] add ForneyLab
```
If you want to be able to use the graph visualization functions, you will also need to have [GraphViz](http://www.graphviz.org/) installed. On Linux, just use `apt-get install graphviz` or `yum install graphviz`. On Windows, run the installer and afterwards manually add the path of the GraphViz installation to the `PATH` system variable. On MacOS, use for example `brew install graphviz`. The `dot` command should work from the command line.

Some demos use the [PyPlot](https://github.com/stevengj/PyPlot.jl) plotting module. Install it using `] add PyPlot`.

Optionally, use `] test ForneyLab` to validate the installation by running the test suite.


## Introduction to ForneyLab

These demos assume that the user is familiar with the FFG formalism. We recommend the following introductions:

1. H.-A. Loeliger, J. Dauwels, J. Hu, S. Korl, Li Ping, and F. Kschischang,
[The factor graph approach to model-based signal processing](https://people.ee.ethz.ch/~papers/docu/aloe-jdau-juhu-skor-2007-1.pdf), Proceedings of the IEEE, vol. 95, no. 6, pp. 1295-1322, June 2007.
2. Korl, Sascha, [A factor graph approach to signal modelling, system identification and filtering](https://www.research-collection.ethz.ch/handle/20.500.11850/82737), Series in signal and information processing
Doctoral Thesis, 2005

We designed ForneyLab to be practical, while retaining maximal flexibility. The inherent modularity of the FFG framework allowed us to make ForneyLab extensible at all levels (nodes, update rules, algorithms, inference engines). Although we had performance in mind while developing ForneyLab, optimally efficient execution of the resulting inference programs (specified in Julia as message passing sequence) may still require custom work.  

The ForneyLab approach to solving inference problems consists of three phases:

1. **Model specification**. ForneyLab provides a simple meta-language to specifiy models.
2. **Message Passing Agorithm (MPA) Generation**. This task is automatically performed by ForneyLab.
3. **MPA Execution**. This is simply evaluating a Julia program.

Each of the demos will step through these phases in turn, showcasing the most important ForneyLab functionalities. For more detailed information we refer to the Julia help functionality (simply type `?` and the ForneyLab function you're interested in), or the source code itself.

## Getting started

There are [demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo) available to get you started.
Additionally, the [ForneyLab project page](http://forneylab.org) contains a talk and other resources that might be helpful.


## License

(c) 2018 GN Store Nord A/S. Permission to use this software for any non-commercial purpose is granted. See `LICENSE.md` file for details.
