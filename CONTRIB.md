Contribution guidelines
=======================


ForneyLab.jl is a probabilistic programming toolbox for message passing-based approximate Bayesian inference. We welcome contributions to ForneyLab.jl from everywhere. This document contains a few development guidelines if you consider contributing.

## Navigating the repository

We organize the ForneyLab repository in the following way:

- `/demo/`: demos in Jupyter (iJulia) notebook format (`.ipynb`)
- `/docs/`: documentation source files
- `/src/`: all source files
    + `algorithms/`: inference algorithm implementations
        * `expectation_propagation/`: EP algorithm implementation
        * `sum_product/`: SP algorithm implementation
        * `variational_bayes/`: VMP algorithm implementation
    + `engines/`: rendering of message passing schedules to executable code
        * `julia/`: Julia engine and update rule implementations
    + `factor_nodes/`: all node-specific files
    + `update_rules/`: message passing update rules
- `/test/`: test files with directory structure similar to `/src/`.

File and directory names always use `snake_case`, except for `REQUIRE` and markdown files in the root directory.

## Reporting bugs

We track bugs using [Github issues](https://github.com/biaslab/ForneyLab.jl/issues). We encourage you to write complete, specific, reproducible bug reports. Mention the versions of Julia and ForneyLab.jl for which you observe unexpected behavior. Please provide a concise description of the problem and complement it with code snippets, test cases, screenshots, tracebacks or any other information that you consider relevant. This will help us to replicate the problem and narrow the search space for solutions.

## Suggesting features

We welcome new feature proposals. However, before submitting a feature request, consider a few things:

- Does the feature require changes in the core ForneyLab.jl code? If it doesn't (for example, you would like to add a node for a particular application), consider making a separate repository for your extensions.
- If you would like to add an implementation of a feature that changes a lot in the core ForneyLab.jl code, please open an issue on Github and describe your proposal first. This will allow us to discuss your proposal with you before you invest your time in implementing something that may be difficult to merge later on.

## Contributing code

### Installing ForneyLab

We suggest that you use the `dev` command from the new Julia package manager to
install ForneyLab.jl for development purposes. To work on your fork of ForneyLab.jl, use your fork's URL address in the `dev` command, for example:

```jl
] dev git@github.com:your_username/ForneyLab.jl.git
```

The `dev` command clones ForneyLab.jl to `~/.julia/dev/ForneyLab`. All local
changes to ForneyLab code will be reflected in imported code.

### Committing code

We use the standard [GitHub Flow](https://guides.github.com/introduction/flow/) workflow where all contributions are added through pull requests. In order to contribute, first [fork](https://guides.github.com/activities/forking/) the repository, then commit your contributions to your fork, and then create a pull request on the `master` branch of the ForneyLab.jl repository.

Before opening a pull request, please make sure that all tests pass without
failing! All demos (can be found in `/demo/` directory) have to run without errors as well.

### Unit tests

We use the test-driven design (TDD) methodology for ForneyLab.jl development. The test coverage should be as complete as possible. Please make sure that you write tests for each piece of code that you want to add.

All unit tests are located in the `/test/` directory. The `/test/` directory follows the structure of the `/src/` directory. Each test file should have following filename format: `test_*.jl`.

The tests can be evaluated by running following command in the Julia REPL:

```jl
] test ForneyLab
```

## Additional materials

ForneyLab.jl is written in the [Julia programming language](julialang.org). In case you are unfamiliar with Julia or would like to deepen your knowledge of Julia, we suggest following resources:

- [Official Julia documentation](https://docs.julialang.org/en/v1/)
- [Learn Julia in Y minutes](https://learnxinyminutes.com/docs/julia/)
- [Think Julia: How to Think Like a Computer Scientist](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html)

If you would like to learn about Git, read for instance the
[ProGit](https://git-scm.com/book/en/v2) book.

We also advise you to use the [Revise.jl](https://github.com/timholy/Revise.jl) package since it reduces the need to restart the Julia REPL whenever you make changes to ForneyLab.jl code.
