ForneyLab.jl
============

Forney-style Factor Graph framework in Julia.

Installation
============

Add the ForneyLab package to your Julia installation:
```jl
Pkg.clone("https://github.com/spsbrats/ForneyLab.jl.git")
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
Once imported, one can instantiate nodes and edges to build a factor graph. It is also possible to create custom nodes and distribution types. When the factor graph has been defined, message passing can be performed on the graph. 

There are [demos](https://github.com/spsbrats/ForneyLab.jl/wiki/ForneyLab-demos) available to get you started. 
More details can be found in the [documentation](http://spsbrats.github.io/ForneyLab/documentation/).
