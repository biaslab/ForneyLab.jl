export huge, tiny, cholinv, diageye, format, *, .*, ^

import Base: *, .*, ^, ==, sqrt

"""ensureMatrix: cast input to a Matrix if necessary"""
ensureMatrix{T<:Number}(arr::AbstractMatrix{T}) = arr
ensureMatrix{T<:Number}(arr::Vector{T}) = Diagonal(arr)
ensureMatrix(n::Number) = fill!(Array(typeof(n),1,1), n)
ensureMatrix(n::Void) = nothing

# Constants to define smallest/largest supported numbers.
# Used for clipping quantities to ensure numerical stability.
const huge = 1e12
const tiny = 1e-12

# Operations related to diagonal matrices
cholinv(m::Number) = 1.0/m
cholinv(M::Matrix) = inv(cholfact(Hermitian(M)))
cholinv(D::Diagonal) = Diagonal(1./D.diag)
diageye(dims::Int64) = Diagonal(ones(dims))

.*(D1::Diagonal, D2::Diagonal) = Diagonal(D1.diag.*D2.diag)
.*(D1::Matrix, D2::Diagonal) = Diagonal(diag(D1).*D2.diag)
.*(D1::Diagonal, D2::Matrix) = D2.*D1

^(D::Diagonal, p::Float64) = Diagonal(D.diag.^p)
sqrt(D::Diagonal) = Diagonal(sqrt(D.diag))

# Symbol concatenation
*(sym::Symbol, num::Number) = Symbol(string(sym, num))
*(num::Number, sym::Symbol) = Symbol(string(num, sym))
*(sym1::Symbol, sym2::Symbol) = Symbol(string(sym1, sym2))
*(sym::Symbol, rng::Range) = Symbol[sym*k for k in rng]
*(rng::Range, sym::Symbol) = Symbol[k*sym for k in rng]
*{T<:Number}(sym::Symbol, vec::Vector{T}) = Symbol[sym*k for k in vec]
*{T<:Number}(vec::Vector{T}, sym::Symbol) = Symbol[k*sym for k in vec]
*(sym1::Symbol, sym2::Vector{Symbol}) = Symbol[sym1*s2 for s2 in sym2]
*(sym1::Vector{Symbol}, sym2::Symbol) = Symbol[s1*sym2 for s1 in sym1]

format(d::Bool) = string(d)

format(d::Symbol) = string(d)

function format(d::Float64)
    if 0.01 < d < 100.0 || -100 < d < -0.01 || d==0.0
        return @sprintf("%.2f", d)
    else
        return @sprintf("%.2e", d)
    end
end

function format{T<:Union{Float64, Bool}}(d::Vector{T})
    s = "["
    for d_k in d[1:end-1]
        s*=format(d_k)
        s*=", "
    end
    s*=format(d[end])
    s*="]"
    return s
end

function format(d::Matrix{Float64})
    s = "["
    for r in 1:size(d, 1)
        s *= format(vec(d[r,:]))
    end
    s *= "]"
    return s
end

format(d::Diagonal{Float64}) = "diag$(format(d.diag))"

function format(v::Vector{Any})
    str = ""
    for (i, entry) in enumerate(v)
        name = replace("$(entry)", "ForneyLab.", "")
        if i < length(v)
            str *= "`$(name)`, "
        else
            str *= "and `$(name)`."
        end
    end
    return str
end

"""isApproxEqual: check approximate equality"""
isApproxEqual(arg1, arg2) = maximum(abs(arg1-arg2)) < tiny

"""isRoundedPosDef: is input matrix positive definite? Round to prevent fp precision problems that isposdef() suffers from."""
isRoundedPosDef(arr::AbstractMatrix{Float64}) = ishermitian(round(Matrix(arr), round(Int, log10(huge)))) && isposdef(Matrix(arr), :L)

function viewFile(filename::AbstractString)
    # Open a file with the application associated with the file type
    is_windows() ? run(`cmd /c start $filename`) : (is_apple() ? run(`open $filename`) : (is_linux() ? run(`xdg-open $filename`) : error("Cannot find an application for $filename")))
end

"""
trigammaInverse(x): solve `trigamma(y) = x` for `y`.

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
Duplicate a method definition with the order of the first two arguments swapped.
This macro is used to duplicate methods that are symmetrical in their first two input arguments,
but require explicit definitions for the different argument orders.
Example:

    @symmetrical function prod!(x, y, z)
        ...
    end
"""
macro symmetrical(orig::Expr)
    (orig.args[1].head == :call) || error("Invalid use of @symmetrical")
    eval(orig)
    swap_arg_indexes = Int64[]
    for i=1:length(orig.args[1].args)
        (typeof(orig.args[1].args[i]) == Expr) || continue
        (orig.args[1].args[i].head == :(::)) || continue
        push!(swap_arg_indexes, i)
        (length(swap_arg_indexes) < 2) || break
    end

    mirrored = deepcopy(orig)
    mirrored.args[1].args[swap_arg_indexes] = orig.args[1].args[reverse(swap_arg_indexes)]
    eval(mirrored)
end

"""
`leaftypes(datatype)` returns all subtypes of `datatype` that are leafs in the type tree.
"""
function leaftypes(datatype::DataType)
    leafs = DataType[]
    stack = [datatype]
    while !isempty(stack)
        for T in subtypes(pop!(stack))
            if isleaftype(T)
                push!(leafs, T)
            else
                push!(stack, T)
            end
        end
    end

    return leafs
end