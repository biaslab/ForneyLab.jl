# Style guide and coding principles for ForneyLab.jl

This document describes style and coding conventions that we use when developing
ForneyLab.jl. We strongly believe that having consistent code style facilitates
contributions and makes code review process easier for all parties.

This style guide is largely based on default [Julia style
guidelines](https://docs.julialang.org/en/v1/manual/style-guide/index.html) and
[PEP-8](https://www.python.org/dev/peps/pep-0008/). Unfortunately, tools like
[rustfmt](https://github.com/rust-lang/rustfmt) and [gofmt](https://golang.org/cmd/gofmt/) are not yet available for Julia.

At the moment not all ForneyLab.jl code follows the rules described in this
guide. We kindly ask to fix any style guide violations in the code that
surrounds your contributions to ForneyLab.jl. If the codebase needs significant
changes sue to style inconsistencies, consider opening a separate pull request
for these changes.

## Style guide

### Line length

We do not set a strict upper limit on number of characters in a line of code. We
believe that imposing a strict limit ends up being counterproductive in cases
where a slightly longer line is more readable then the alternatives, but is
reformatted due to strict coding conventions.

Instead we suggest following approach to organizing your code. Imagine three
zones: green, yellow, and red. Green is under 80 characters, yellow is between
80 and 100 characters and red is above 100 characters. Try to write your code
in such a way that ~95% of lines of your code are in the green zone, and no code
is in the red zone.

### Method definitions

If a method definition doesn't fit in a single line of code, always use
multiline form of method definition:

```julia
# Yes:
function fn(item::Float64) = sqrt(item) + 2.71

# No:
function fn(data::Vector{T}, item::T) where {T<:Float64} = [sqrt(item) + exp(abs(x)) + 2.71 for x in data] 

# Yes:
function fn(data::Vector{T}, item::T) where T<:Float64
    return [sqrt(item) + exp(abs(x)) + 2.71 for x in data] 
end
```

If method definition (method name + parameter line) exceeds character limit,
separate each parameter by a newline and indent once:

```julia
# Yes: 
function condense(schedule::Schedule)
    ...
end

# No: 
function condense(schedule::Schedule, fg::FactorGraph, variable::Symbol, id=generateId(variable))
    ...
end

# Yes:
function condense(
    schedule::Schedule, 
    fg::FactorGraph, 
    variable::Symbol, 
    id=generateId(variable)
)
    ...
end

```


### `return` keyword

Always use `return` keyword in multiline function definitions:

```julia
# Yes:
function fn(a::Int)
    return a^2
end

# No:
function fn(a::Int)
    a^2
end
```

### `using` and `import`

Whenever possible, import only particular things in the scope. E.g., rather than
`using SpecialFunctions`, write `using SpecialFunctions: lgamma`. This way you
avoid bringing unnecessary objects into scope.

Use `import` keyword only if you extend the imported methods: `import Base: show`.