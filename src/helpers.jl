export huge, tiny, cholinv, diageye, eye, format, *, ^, mat, step!, init

"""Cast input to a `Matrix` if necessary"""
ensureMatrix(arr::AbstractMatrix{T}) where T<:Number = arr
ensureMatrix(arr::Vector{T}) where T<:Number = Diagonal(arr)
ensureMatrix(n::Number) = fill!(Array{typeof(n)}(undef,1,1), n)
ensureMatrix(n::Nothing) = nothing

# Constants to define smallest/largest supported numbers.
# Used for clipping quantities to ensure numerical stability.
const huge = 1e12
const tiny = 1e-12

"""
Wrapper for `logabsgamma` function that returns first element of its output
"""
function labsgamma(x::Number)
    return logabsgamma(x)[1]
end

"""
Wrapper for `logabsbeta` function that returns first element of its output
"""
function labsbeta(x::Number, y::Number)
    return logabsbeta(x, y)[1]
end


"""
Matrix inversion using Cholesky decomposition,
attempts with added regularization (1e-8*I) on failure.
"""
function cholinv(M::AbstractMatrix)
    try
        return inv(cholesky(Hermitian(Matrix(M))))
    catch
        try
            return inv(cholesky(Hermitian(Matrix(M) + 1e-8*I)))
        catch exception
            if isa(exception, PosDefException)
                error("PosDefException: Matrix is not positive-definite, even after regularization. $(typeof(M)):\n$M")
            else
                println("cholinv() errored when inverting $(typeof(M)):\n$M")
                rethrow(exception)
            end
        end
    end
end
cholinv(m::Number) = 1.0/m
cholinv(D::Diagonal) = Diagonal(1 ./ D.diag)
eye(n::Number) = Diagonal(I,n)
diageye(dims::Int64) = Diagonal(ones(dims))

# Symbol concatenation
*(sym::Symbol, num::Number) = Symbol(string(sym, num))
*(num::Number, sym::Symbol) = Symbol(string(num, sym))
*(sym1::Symbol, sym2::Symbol) = Symbol(string(sym1, sym2))
*(sym::Symbol, rng::AbstractRange) = Symbol[sym*k for k in rng]
*(rng::AbstractRange, sym::Symbol) = Symbol[k*sym for k in rng]
*(sym::Symbol, vec::Vector{T}) where T<:Number = Symbol[sym*k for k in vec]
*(vec::Vector{T}, sym::Symbol) where T<:Number = Symbol[k*sym for k in vec]
*(sym1::Symbol, sym2::Vector{Symbol}) = Symbol[sym1*s2 for s2 in sym2]
*(sym1::Vector{Symbol}, sym2::Symbol) = Symbol[s1*sym2 for s1 in sym1]

format(d::Bool) = string(d)

format(d::Symbol) = string(d)

format(d::Int) = string(d)

function format(d::AbstractFloat)
    if 0.01 < d < 100.0 || -100 < d < -0.01 || d==0.0
        return @sprintf("%.2f", d)
    else
        return @sprintf("%.2e", d)
    end
end

function format(d::Vector{T}) where T<:Union{AbstractFloat, Bool}
    s = "["
    for d_k in d[1:end-1]
        s*=format(d_k)
        s*=", "
    end
    s*=format(d[end])
    s*="]"
    return s
end

function format(d::Matrix{T}) where T<:AbstractFloat
    s = "["
    for r in 1:size(d, 1)
        s *= format(vec(d[r,:]))
    end
    s *= "]"
    return s
end

format(d::Diagonal{T}) where T<:AbstractFloat = "diag$(format(d.diag))"

function format(v::Vector{Any})
    str = ""
    for (i, entry) in enumerate(v)
        name = replace("$(entry)", "ForneyLab." => "")
        if i < length(v)
            str *= "`$(name)`, "
        else
            str *= "and `$(name)`."
        end
    end
    return str
end

"""Check if arguments are approximately equal"""
isApproxEqual(arg1, arg2) = maximum(abs.(arg1-arg2)) < tiny

"""Checks if input matrix is positive definite. We also perform rounding in order to prevent floating point precision problems that `isposdef()`` suffers from."""
isRoundedPosDef(arr::AbstractMatrix{Float64}) = ishermitian(round.(Matrix(arr), digits=round(Int, log10(huge)))) && isposdef(Hermitian(Matrix(arr), :L))

function viewFile(filename::AbstractString)
    # Open a file with the application associated with the file type
    Sys.iswindows() ? run(`cmd /c start $filename`) : (Sys.isapple() ? run(`open $filename`) : (is_linux() ? run(`xdg-open $filename`) : error("Cannot find an application for $filename")))
end

"""
Solve `trigamma(y) = x` for `y`.

Uses Newton's method on the convex function 1/trigramma(y).
Iterations converge monotonically.
Based on trigammaInverse implementation in R package "limma" by Gordon Smyth:
https://github.com/Bioconductor-mirror/limma/blob/master/R/fitFDist.R
"""
function trigammaInverse(x::Float64)
    y = 0.5 + 1/x
    iter = 0
    while iter < 20
        iter += 1
        tri = trigamma(y)
        dif = tri*(1-tri/x) / polygamma(2, y)
        y += dif
        ((-dif/y) > 1e-8) || break
    end

    return y
end

"""
    @symmetrical `function_definition`

Duplicate a method definition with the order of the first two arguments swapped.
This macro is used to duplicate methods that are symmetrical in their first two input arguments,
but require explicit definitions for the different argument orders.
Example:

    @symmetrical function prod!(x, y, z)
        ...
    end
"""
macro symmetrical(fn::Expr)
    # Check if macro is applied to a function definition
    # Valid function definitions include:
    # 1. foo([ args... ]) [ where ... [ where ... [ ... ] ] ] = :block
    # 2. function foo([ args... ]) [ where ... [ where ... [ ... ] ] ]
    #        :block
    #    end
    if (fn.head === :(=) || fn.head === :function) &&
        (fn.args[1] isa Expr && fn.args[2] isa Expr) &&
        (fn.args[2].head === :block)
        return esc(quote
            $fn
            $(swap_arguments(fn))
        end)
    else
        error("@symmetrical macro can be applied only to function definitions")
    end
end

function swap_arguments(fn::Expr)
    swapped = copy(fn)

    if swapped.args[1].head === :where
        swapped.args[1] = swap_arguments(swapped.args[1])
    elseif swapped.args[1].head === :call && length(fn.args[1].args) >= 3 # Note: >= 3, because the first argument is a function name
        swapped.args[1].args[2] = fn.args[1].args[3]
        swapped.args[1].args[3] = fn.args[1].args[2]
    else
        error("Function method passed for @symmetrical macro must have more than 2 arguments")
    end

    return swapped
end

"""
`leaftypes(datatype)` returns all subtypes of `datatype` that are leafs in the type tree.
"""
function leaftypes(datatype::Type)
    leafs = []
    stack = Type[datatype]
    # push!(stack, datatype)
    while !isempty(stack)
        for T in subtypes(pop!(stack))
            if isconcretetype(T)
                push!(leafs, T)
            else
                push!(stack, T)
            end
        end
    end

    return leafs
end

"""
Helper function to construct 1x1 `Matrix`
"""
mat(sc) = reshape([sc],1,1)

"""
Helper function to call dynamically generated `step!` functions
"""
function step!(id::Symbol, args...)
    getfield(Main, :step*id*:!)(args...)
end

"""
Helper function to call dynamically generated `init` functions
"""
function init(id::Symbol)
    getfield(Main, :init*id)()
end
