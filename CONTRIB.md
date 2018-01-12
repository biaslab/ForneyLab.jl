Contribution guidelines
=======================

File structure
--------------
- `/demo/`: demos in Jupyter (iJulia) notebook format (`.ipynb`)
- `/src/`: all source files
    + `algorithms/`: algorithm implementations
        * `expectation_propagation/`: EP algorithm implementation
        * `sum_product/`: SP algorithm implementation
        * `variational_bayes/`: VMP algorithm implementation
    + `engines/`: render message passing schedules to executable code
        * `julia/`: Julia engine and update rule implementations
    + `factor_nodes/`: all node-specific files
    + `update_rules/`: message passing update rules 
- `/test/`: test files with similar directory structure as `/src/`.

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
    mv $(basename "$s") $s
    git add $s
done
```


Style conventions
-----------------
We use the default [Julia style conventions](http://julia.readthedocs.org/en/latest/manual/style-guide/), with a few modifications. The most important points (and modifications to the Julia style) are:

- Indentation: 4 spaces
- Type names in `UpperCamelCase`
- Function names in `lowerCamelCase` (differs from the official Julia convention)
- Variable names and function arguments in `snake_case`
- File and directory names in `snake_case`
- The name of a method that writes back in its argument(s) end in `!`
- Method headers (comments) come after the method definition and are indented together with the method implementation if the method implementation takes more than 1 line.


Testing setup
-------------
The `/test/` directory follows the structure of the `/src/` directory. Every test file should have the following filename format: `test_*.jl`.
The test coverage should be as complete as possible. Please make sure you write tests for every piece of code you add.

**Before each commit, confirm that all tests pass**. Pull requests that break tests or that do not include sufficient tests will not be merged. The tests are evaluated by simply running

```jl
Pkg.test("ForneyLab")
```