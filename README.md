ForneyLab.jl
============

ForneyLab.jl is a Julia package for automatic generation of (Bayesian) inference algorithms. Given a probabilistic model, ForneyLab generates efficient Julia code for message-passing based inference. It uses the model structure to generate an algorithm that consists of a sequence of local computations on a Forney-style factor graph (FFG) representation of the model. For an excellent introduction to message passing and FFGs, see [The Factor Graph Approach to Model-Based Signal Processing](https://ieeexplore.ieee.org/document/4282128/) by Loeliger et al.

We designed ForneyLab with a focus on flexible and modular modeling of time-series data. ForneyLab enables a user to:

- Conveniently specify a probabilistic model;
- Automatically generate an efficient inference algorithm;
- Compile the inference algorithm to executable Julia code.

The current version supports [belief propagation](https://en.wikipedia.org/wiki/Belief_propagation) (sum-product message passing), [variational message passing](https://en.wikipedia.org/wiki/Variational_message_passing) and [expectation propagation](https://en.wikipedia.org/wiki/Expectation_propagation).

For a more detailed introduction we refer to the [demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo).

Installation
============

Add the ForneyLab package to your Julia installation:
```jl
Pkg.clone("https://github.com/biaslab/ForneyLab.jl.git")
```
To update, use:
```jl
Pkg.update()
```
If you want to be able to use the graph visualization functions, you'll also need to have [GraphViz](http://www.graphviz.org/) installed. On Linux, just use `apt-get install graphviz` or `yum install graphviz`. On Windows, run the installer and afterwards manually add the path of the GraphViz installation to the `PATH` system variable. The command `dot -?` should work from the command line.

Some demos use the [PyPlot](https://github.com/stevengj/PyPlot.jl) plotting module. Install it using `Pkg.add("PyPlot")`.

Getting started
===============

Import ForneyLab:
```jl
using ForneyLab
```

There are [demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo) available to get you started.

License
=======

(c) 2018 GN Store Nord A/S. Permission to use this software for any non-commercial purpose is granted. See `LICENSE.md` file for details.