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
To update an already installed package, use:
```jl
Pkg.update()
```
Usage
=====
Import ForneyLab:
```jl
using ForneyLab
```
Once imported, one can create nodes and edges to build a factor graph, or create custom node and message types. There are demos available [online](http://192.71.151.86/ForneyLab.jl-demos/) to get you started.

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