ForneyLab.jl
============

Forney-style Factor Graph framework in Julia.
**This software is still in an early development stage**.

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

Usage
=====
Import ForneyLab:
```jl
using ForneyLab
```
Once imported, one can create nodes and edges to build a factor graph, or create custom node and message types.

Getting started
===============

There are [demos](https://github.com/spsbrats/ForneyLab.jl/wiki/ForneyLab-demos) available to get you started (which can also be accessed directly through [here](http://192.71.151.86/ForneyLab.jl-demos/)).

Extending ForneyLab.jl
======================

It is easy to extend ForneyLab with for example custom node types and message types.
Extensions live in a separate directory, and are implemented in the following way:

- Create a directory to hold your extensions. This directory can be anywhere, and can be a separate git repo.
- The extension directory should at least contain the directory `src` and the file `src/forneylab_extensions.jl`. Ideally, one would use the same directory tree as the ForneyLab package itself. So, a custom node would be implemented in a separate file `[extension-dir]/src/nodes/my_node.jl` and this file should be included in `[extension-dir]/src/forneylab_extensions.jl` by adding `include("nodes/my_node.jl")` to that file. It is recommended to also add a main test file `[extension-dir]/test/test_forneylab_extensions.jl` for your extensions. This way, your extensions will be tested together with ForneyLab itself.
- To use your extensions, define `FORNEYLAB_EXTENSION_DIR` in your top-level script before including ForneyLab itself:

```jl
FORNEYLAB_EXTENSION_DIR = "/home/marco/dev/my_extensions"
using ForneyLab
```
- When you run ForneyLab's tests, the tests in `FORNEYLAB_EXTENSION_DIR/test/test_forneylab_extensions.jl` will automatically be evaluated (if this file exists). Moreover, the general ForneyLab tests will also be evaluated on your custom extensions (i.e. node types).

Troubleshooting
---------------
When you run into trouble installing IJulia and the package manager begins to complain about not being able to build LibCURL and Nettle, first take a look at [this gist](https://gist.github.com/ThijsvdLaar/8e8f48077e5373ab7b80).
