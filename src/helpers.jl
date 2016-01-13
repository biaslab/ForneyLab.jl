export isApproxEqual, huge, tiny, rules, format, isValid, invalidate!, *

import Base.*, Base.==
# ensureMatrix: ensure that the input is a 2D array or nothing
ensureMatrix{T<:Number}(arr::Array{T, 2}) = arr
ensureMatrix{T<:Number}(arr::Array{T, 1}) = reshape(arr, 1, 1)
ensureMatrix(n::Number) = fill!(Array(typeof(n),1,1), n)
ensureMatrix(n::Void) = nothing

# Constants to define smallest/largest supported numbers.
# Used for clipping quantities to ensure numerical stability.
const huge = 1e12
const tiny = 1e-12

# Functions for checking validity and invalidating arrays of floats
# An array is invalid if and only if its first entry is NaN
isValid(v::Array{Float64}) = !isnan(v[1])
isValid(v::Float64) = !isnan(v)

function invalidate!(v::Array{Float64})
    v[1] = NaN
    return v
end

# Symbol concatenation
*(sym::Symbol, num::Number) = symbol(string(sym, num))
*(num::Number, sym::Symbol) = symbol(string(num, sym))
*(sym1::Symbol, sym2::Symbol) = symbol(string(sym1, sym2))

function format(d::Bool)
    string(d)
end

function format(d::Float64)
    if 0.01 < d < 100.0 || -100 < d < -0.01 || d==0.0
        return @sprintf("%.2f", d)
    else
        return @sprintf("%.2e", d)
    end
end

function format(d::Vector{Float64})
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
    for r in 1:size(d)[1]
        s *= format(vec(d[r,:]))
    end
    s *= "]"
end

# isApproxEqual: check approximate equality
isApproxEqual(arg1, arg2) = maximum(abs(arg1-arg2)) < tiny

# isRoundedPosDef: is input matrix positive definite? Round to prevent fp precision problems that isposdef() suffers from.
isRoundedPosDef{T<:AbstractFloat}(arr::Array{T, 2}) = ishermitian(round(arr, round(Int, log(10, huge)))) && isposdef(arr, :L)

function viewFile(filename::AbstractString)
    # Open a file with the application associated with the file type
    @windows? run(`cmd /c start $filename`) : (@osx? run(`open $filename`) : (@linux? run(`xdg-open $filename`) : error("Cannot find an application for $filename")))
end

function truncate(str::ASCIIString, max::Integer)
    # Truncate srt to max positions
    if length(str)>max
        return "$(str[1:max-3])..."
    end
    return str
end

function pad(str::ASCIIString, size::Integer)
    # Pads str with spaces until its length reaches size
    str_trunc = truncate(str, size)
    return "$(str_trunc)$(repeat(" ",size-length(str_trunc)))"
end

immutable HTMLString
    s::ByteString
end
import Base.writemime
writemime(io::IO, ::MIME"text/html", y::HTMLString) = print(io, y.s)

function rules(node_type::Union{DataType, Void}=nothing; format=:table)
    # Prints a table or list of node update rules
    rule_dict = YAML.load_file("$(Pkg.dir("ForneyLab"))/src/update_equations.yaml")
    all_rules = rule_dict["rules"]

    # Select node specific rules
    if node_type != nothing
        rule_list = Dict()
        for (id, rule) in all_rules
            if rule["node"] == "$(node_type)"
                merge!(rule_list, Dict{Any,Any}(id=>rule)) # Grow the list of rules
            end
        end
    else
        rule_list = all_rules
    end

    node="node"; reference="reference"; formula="formula"; diagram="diagram"
    # Write rule list to output
    if format==:table
        println("Message calculation rules (node id (reference))")
        println("-----------------------------------------------")
        for (id, rule) in rule_list
            println("$(rule[node]) $(id) ($(rule[reference]))")
        end
        println("\nUse rules(NodeType) to view all rules of a specific node type; use rules(..., format=:list) to view the formulas in latex output.")
    elseif format==:list
        display(HTMLString("<p><b>Parameterizations</b><br>Student's-t distribution: St(m,W,ν)<br>Delta Distribution: δ(m)<br>Gaussian distribution (mean, variance): N(m,V), N(m,W⁻¹) for multivariate, or N(m,λ⁻¹) for univariate<br>Gamma distribution (shape, rate): Gam(a,b)<br>Inverse gamma distribution (shape, scale): Ig(a,b)<br>Beta distribution: Bet(a,b)<br>Normal-gamma distribution: Ng(m,β,a,b)</p>"))
        display(HTMLString("<p>Parameters in the equations are denoted by the parameter character (see above list for distribution parametrizations) and a subcript with the variable it applies to.</p>"))
        for (id, rule) in rule_list
            display(HTMLString("<b>$(id)</b>; $(rule[reference]):"))
            display(HTMLString("<font face=\"monospace\">$(rule[diagram])</font>"))
            display(LaTeXString(rule[formula]))
            display(HTMLString("<br>"))
        end
    else
        error("Unknown format $(format), use :table or :list instead")
    end
end
