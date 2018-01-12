ForneyLab.jl
============

Forney-style Factor Graph framework in Julia.

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
