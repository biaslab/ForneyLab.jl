# Style guide and coding principles for ForneyLab.jl

This document describes style and coding conventions that we use when developing
ForneyLab.jl. We strongly believe that having consistent code style facilitates
contributions and makes code review process easier for all parties.

This style guide is largely based on default Julia style guidelines and PEP-8.

Not all ForneyLab.jl code follows the rules described in this guide. We kindly
ask to fix any style guide violations in the code that surrounds your
contributions to ForneyLab.jl. If the codebase needs significant changes,
consider opening a separate pull request for these changes.

## Style guide

### Line length

We do not set a hard upper limit on number of characters in a line of code.
However, we suggest that you 

### Method definitions

If method definition doesn't fit in a single line of code, always use
multiline form of function definition:

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

If method's signature (parameter line) exceeds soft character limit, separate
each parameter by a newline and indent once:

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
