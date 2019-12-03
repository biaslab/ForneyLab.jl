# Style guide and coding principles for ForneyLab.jl

This document describes style and coding conventions that we use when developing
ForneyLab.jl. We strongly believe that having consistent code style facilitates
contributions and makes code review process easier for all parties.
Unfortunately, tools like [rustfmt](https://github.com/rust-lang/rustfmt) and
[gofmt](https://golang.org/cmd/gofmt/) are not yet available for Julia. For this
reason we have assembled a set of guidelines for writing source code and
documentation.

This style guide is largely inspired by [Julia style
guidelines](https://docs.julialang.org/en/v1/manual/style-guide/index.html) and
[PEP-8](https://www.python.org/dev/peps/pep-0008/). It is also by no means
comprehensive. We suggest that readers get familiar with standard Julia style
guidelines first as we do not repeat conventions that are listed there (e.g.
appending `!` to names of the methods that modify their arguments).

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
80 and 100 characters and red is above 100 characters. Try to write your code in
such a way that ~95% of your code is in the green zone, and no code is in the
red zone.

### Indentation

Julia, unlike Python, is largely insensitive to whitespace characters and
indentation levels. However, for consistency, use 4 spaces per indentation
level. Start all new blocks with new level of indentation except for the first
level inside the module:

```julia
# Yes:
module MyModule

function add2(x::Number)
    return x+2
end

end

# No:
module MyModule

    function add2(x::Number)
    return x+2
    end

end
```

### Whitespace

Avoid extraneous whitespace in the following situations:

- Immediately inside parentheses, brackets or braces.

    ```julia
    # Yes:
    spam(ham[1], Dict("eggs"=>2))
    # No:  
    spam( ham[ 1 ], Dict( "eggs"=>2 ) )
    ```

- Between a trailing comma and a following close parenthesis.
    
    ```julia
    # Yes: 
        foo = (0,)
    # No:  
        bar = (0, )
    ```

- Immediately before a comma, semicolon, or colon:
    
    ```julia
    Yes: if x == 4 println((x, y)); x, y = y, x end
    No:  if x == 4 println((x , y)) ; x , y = y , x end
    ```

- When using ranges (unless additional parameters are used):

    ```julia
    Yes: ham[1:9], ham[1:9:3], ham[1:9:end]
    No: ham[1: 9], ham[1 :9], ham[1:9 :end]
    
    #Yes:
    ham[lower:upper], ham[lower:step:end]
    ham[lower+offset : upper+offset]
    ham[lower + offset : upper + offset]
    ```

- More than one space around an assignment (or other) operator to align it with
  another.

    ```julia
    #Yes:
    x = 1
    y = 2
    long_variable = 3

    #No:
    x             = 1
    y             = 2
    long_variable = 3
    ```

Other rules:

- Always surround these binary operators with a single space on either side:
  assignment (=), updating operators (+=, -= etc.), comparisons (==, ===, <, >,
  !=,
  <>, <=, >=, etc.), Booleans (and, or, not), infix operators (`is`, `isa`).

- If operators with different priorities are used, consider adding whitespace around the operators with the lowest priority(ies). However, never use more than one space, and always have the same amount of whitespace on both sides of a binary operator.

    ```julia
    # Yes:
    i = i + 1
    submitted += 1
    x = x*2 - 1
    hypot2 = x*x + y*y
    c = (a+b) * (a-b)
    c = (a + b)*(a - b)
    
    # No:
    i=i+1
    submitted +=1
    x = x * 2 - 1
    hypot2 = x * x + y * y
    c = (a + b) * (a - b)
    ```

### Naming conventions
- Type names use `UpperCamelCase`: `FactorNode`, `Gaussian`
- Function names are `lowerCamelCase` (differs from the official Julia
  convention): `isApplicable`, `currentGraph`
- Variable names and function arguments (e.g. `inbound_messages`) use `snake_case`

### Method definitions

Use one-line method definitions only if they fit on a single line:

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
    id=generateId(variable))

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

`using`/`import` statements are always put at the top of the file, just after
any module comments and docstrings, and before module globals and constants.

Whenever possible, import only specific things into the scope. E.g., rather than
`using SpecialFunctions`, write `using SpecialFunctions: lgamma`. This way you
avoid bringing unnecessary objects into scope.

Use `import` keyword only if you extend the imported methods: `import Base:
show`.

### Comments

Whenever you write comments, avoid simply explicitly listing what does the code
do:

```julia
# Yes:
i = i+1     # Correct the index for the dummy message

# No:
i = i+1     # Increase i
```

Comments that contradict the code are worse than no comments. Always make a
priority of keeping the comments up-to-date when the code changes!

Comments should be complete sentences. The first word should be capitalized,
unless it is an identifier that begins with a lower case letter (never alter the
case of identifiers!).

Block comments generally consist of one or more paragraphs built out of complete
sentences, with each sentence ending in a period.

### Documentation

All exported types and functions have to have docstrings. We also recommend
including docstrings for most functions and types.

We provide useful documentation template below:

#### Documentation template for a node type

```julia
"""
Description:

    A gamma node with shape-rate parameterization:

    f(out,a,b) = Gam(out|a,b) = 1/Γ(a) b^a out^{a - 1} exp(-b out)

Interfaces:

    1. out
    2. a (shape)
    3. b (rate)

Construction:

    Gamma(out, a, b, id=:some_id)
"""
mutable struct Gamma <: SoftFactor
    ...
end
```

We suggest to include following in the docstring:

- Brief description (1-2 sentences)
- Node function
- List of interfaces and corresponding variables

For more details on documenting Julia code we refer the reader to Julia
[documentation](https://docs.julialang.org/en/v1.1/manual/documentation/).
