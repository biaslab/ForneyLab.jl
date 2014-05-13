Contribution guidelines
=======================

File structure
--------------
- `/demo/`: demos in iJulia notebook format (`.ipynb`)
- `/src/`: all source files
    + `ForneyLab.jl`: top-level module with general definitions and functions
    + `helpers.jl`: generic helper functions functions
    + `messages.jl`: all message type definitions and message-specific functions
    + `nodes/`: contains all node-specific files
        * `composite/`: composite node files
            - `[node_name].jl`: node-specific definitions and methods
        * `[node_name].jl`: node-specific definitions and methods
- `/test/`: FactCheck test files for every file in `/src/` (identical directory tree). File format: `test_[src-filename].jl`

File and directory names are always in `snake_case`, except for `REQUIRE` and markdown files in the root directory.

Style conventions
-----------------
We use the default [Julia style conventions](http://julia.readthedocs.org/en/latest/manual/style-guide/), with a few modifications. The most important points (and modifications to the Julia style) are:

- Indentation: 4 spaces
- Type names in `UpperCamelCase`
- Method names in `lowerCamelCase` (different from Julia default)
- Variable names and function arguments in `snake_case`
- File and directory names in `snake_case`
- The name of a method that writes back in its argument(s) end in `!`
- Method headers (comments) come after the method definition and are indented together with the method implementation if the method implementation takes more than 1 line.

Apart from this, there are some project-specific conventions:

- The name of a subtype of `Node` always ends in `Node`. Example: `EqualityNode`.
- The type name of a composite node always ends in `CompositeNode`. Example: `GainEqualityCompositeNode`.
- The name of a subtype of `Message` always ends in `Message`. Example: `GaussianMessage`.

Testing setup
-------------
[FactCheck](https://github.com/zachallaun/FactCheck.jl) is used as testing framework.
For every `.jl` file in `/src/`, there should be a corresponding FactCheck test file in `/test/`. File format: `test_[src-filename].jl`. The test file should test the behavior of the implementation in the source file as extensively as possible. In general, at least the following should be tested:

- The correctness of all defined constructors
- Correctness of the output of all functions for all valid input types
- Exceptions for invalid input arguments
- Triggering all possible exceptions

The main testing file `test_forneylab.jl` performs some tests on all available node types.
**Before each commit, confirm that all tests pass**. Pull requests that break tests or that do not include sufficient tests will not be merged. The tests are evaluated by simply including the main test file:

```jl
using ForneyLab
include("test/test_forneylab.jl")
```

If you add a new test file, do not forget to add the inclusion to `/test/test/test_forneylab.jl`.
