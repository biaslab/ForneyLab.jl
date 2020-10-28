ForneyLab.jl
============

[![Build Status](https://travis-ci.org/biaslab/ForneyLab.jl.svg?branch=master)](https://travis-ci.org/biaslab/ForneyLab.jl)
[![Documentation](https://img.shields.io/badge/doc-stable-blue.svg)](https://biaslab.github.io/ForneyLab.jl/stable/)

ForneyLab.jl is a Julia package for automatic generation of (Bayesian) inference algorithms. Given a probabilistic model, ForneyLab generates efficient Julia code for message-passing based inference. It uses the model structure to generate an algorithm that consists of a sequence of local computations on a Forney-style factor graph (FFG) representation of the model. For an excellent introduction to message passing and FFGs, see [The Factor Graph Approach to Model-Based Signal Processing](https://ieeexplore.ieee.org/document/4282128/) by Loeliger et al. (2007). Moreover, for a comprehensive overview of the underlying principles behind this tool, see [A Factor Graph Approach to Automated Design of Bayesian Signal Processing Algorithms](https://arxiv.org/abs/1811.03407) by Cox et. al. (2018).

We designed ForneyLab with a focus on flexible and modular modeling of time-series data. ForneyLab enables a user to:

- Conveniently specify a probabilistic model;
- Automatically generate an efficient inference algorithm;
- Compile the inference algorithm to executable Julia code.

The current version supports [belief propagation](https://en.wikipedia.org/wiki/Belief_propagation) (sum-product message passing), [variational message passing](https://en.wikipedia.org/wiki/Variational_message_passing) and [expectation propagation](https://en.wikipedia.org/wiki/Expectation_propagation).

The [ForneyLab project page](http://forneylab.org) provides more background on ForneyLab as well as pointers to related literature and talks. For a practical introduction, have a look at [the demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo).

Documentation
=============

Full documentation is available at [BIASlab website](https://biaslab.github.io/forneylab/docs/).

It is also possible to build documentation locally. Just execute

```bash
$ julia make.jl
```

in the `docs/` directory to build a local version of the documentation.

Installation
============

Install ForneyLab through the Julia package manager:
```jl
] add ForneyLab
```
If you want to be able to use the graph visualization functions, you will also need to have [GraphViz](http://www.graphviz.org/) installed. On Linux, just use `apt-get install graphviz` or `yum install graphviz`. On Windows, run the installer and afterwards manually add the path of the GraphViz installation to the `PATH` system variable. On MacOS, use for example `brew install graphviz`. The `dot` command should work from the command line.

Some demos use the [PyPlot](https://github.com/stevengj/PyPlot.jl) plotting module. Install it using `] add PyPlot`.

Optionally, use `] test ForneyLab` to validate the installation by running the test suite.


Getting started
===============

There are [demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo) available to get you started.
Additionally, the [ForneyLab project page](http://forneylab.org) contains a talk and other resources that might be helpful.


License
=======

(c) 2019 GN Store Nord A/S. Permission to use this software for any non-commercial purpose is granted. See `LICENSE.md` file for details.
