export isApproxEqual, huge, tiny, rules

# ensureMatrix: ensure that the input is a 2D array or nothing
ensureMatrix{T<:Number}(arr::Array{T, 2}) = arr
ensureMatrix{T<:Number}(arr::Array{T, 1}) = reshape(arr, 1, 1)
ensureMatrix(n::Nothing) = nothing

# Define what is considered big and small, mostly used for setting vague messages.
huge(::Type{Float64}) = 1e12
huge() = huge(Float64)
tiny(::Type{Float64}) = 1./huge(Float64)
tiny() = tiny(Float64)

# isApproxEqual: check approximate equality
isApproxEqual(arg1, arg2) = maximum(abs(arg1-arg2)) < tiny()

# isRoundedPosDef: is input matrix positive definite? Round to prevent fp precision problems that isposdef() suffers from.
isRoundedPosDef{T<:FloatingPoint}(arr::Array{T, 2}) = ishermitian(round(arr, int(log(10, huge())))) && isposdef(arr, :L)

function viewFile(filename::String)
    # Open a file with the application associated with the file type
    @windows? run(`cmd /c start $filename`) : (@osx? run(`open $filename`) : (@linux? run(`xdg-open $filename`) : error("Cannot find an application for $filename")))
end

global unnamed_counter = 0
function unnamedStr()
    global unnamed_counter::Int64 += 1
    return "unnamed$(unnamed_counter)"
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

type Latex
    s::String
end

import Base.writemime
writemime(io::IO, ::MIME"text/latex", y::Latex) = print(io, y.s)

function rules(node_type::Union(DataType, Nothing)=nothing; format=:table)
    # Prints a table or list of node update rules
    rule_dict = YAML.load(open("Code/ForneyLab.jl/src/update_equations.yaml"))
    all_rules = rule_dict["rules"]
    
    # Select node specific rules
    if node_type != nothing
        rule_list = Dict[]
        for rule in all_rules
            if rule["node"] == "$(node_type)"
                push!(rule_list, rule) # Grow the list of rules
            end
        end
    else
        rule_list = all_rules
    end

    id="id"; reference="reference"; formula="formula"
    # Write rule list to output
    if format==:table
        println("|              id              |                    reference                     |                    formula                    |")
        println("|------------------------------|--------------------------------------------------|-----------------------------------------------|")
        for rule in rule_list
            println("|$(pad(rule[id],30))|$(pad(rule[reference],50))|$(pad(rule[formula],47))|")
        end
        println("\nUse rules(NodeType) to view all rules of a specific node type; use rules(..., format=:list) to view the formulas in latex output.")
    elseif format==:list
        for rule in rule_list
            display(Latex("$(rule[id]); $(rule[reference]):"))
            display(Latex(rule[formula]))
        end
    else
        error("Unknown format $(format), use :table or :list instead")
    end
end