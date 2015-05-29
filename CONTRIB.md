Contribution guidelines
=======================

File structure
--------------
- `/demo/`: demos in Jupyter (iJulia) notebook format (`.ipynb`)
    + `images/`: images required in demo's
        * `src/`: latex source files for images
- `/src/`: all source files
	+ `algorithms/`: algorithm implementations
    + `distributions/`: distribution types and functions, including marginal calculations
    + `nodes/`: contains all node-specific files
        * `composite/`: composite node files
- `/test/`: FactCheck test files with identical directory structure as `/src/`.

File and directory names are always in `snake_case`, except for `REQUIRE` and markdown files in the root directory.

Pre-commit hook for demos
-------------------------

The Jupyter notebook demos should always be committed in an executed state. To make sure that a new or modified notebook is fully executed when it is committed, one might use a pre-commit hook. The script below will execute all staged demo notebooks which have been added or modified. Create file `.git/hooks/pre-commit` in your local repo, and paste the following shell script to enable the pre-commit hook:

```
#!/bin/bash
#
# Execute all staged added/modified Jupyter(iJulia) demos
# Requires iPython v3.0 or higher

for s in `git diff --staged --name-only --diff-filter=AM | grep "demo/.*\.ipynb"`; do 
	echo "Executing added or modified demo: " $s
	ipython nbconvert --to=notebook --execute --inplace $s
done
```


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
- The type name of a composite node always ends in `CompositeNode`. Example: `GainEqualityNode`.
- The name of a subtype of `ProbabilityDistribution` always ends in `Distribution`. Example: `GaussianDistribution`.

Testing setup
-------------
[FactCheck](https://github.com/zachallaun/FactCheck.jl) is used as testing framework. The `/test/` directory follows the structure of the `/src/` directory.

In general, at least the following should be tested:

- The correctness of all defined constructors
- Correctness of the output of all functions for all valid input types
- Exceptions for invalid input arguments
- Triggering all possible exceptions

The main testing file `test_forneylab.jl` performs some tests on all available node types.
**Before each commit, confirm that all tests pass**. Pull requests that break tests or that do not include sufficient tests will not be merged. The tests are evaluated by simply running or including `test/runtests.jl`:

```jl
include("test/runtests.jl")
```
Alternatively, you can also use:

```jl
Pkg.test("ForneyLab")
```

If you add a new test file, do not forget to add the inclusion to `/test/test_forneylab.jl`.

To test the demos, use:

```jl
include("test/test_demos.jl")
```
