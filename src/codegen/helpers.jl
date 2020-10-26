function get_options(options_expr::Expr)
    options = Dict{Symbol,Any}()
    
    options_expr.head == :vect || return :(error("Incorrect options specification: options argument must be a vector expression"))
    
    for arg in options_expr.args
        arg isa Expr && arg.head == :(=) || return :(error("Incorrect options specification: options item must be an assignment expression"))
        options[arg.args[1]] = arg.args[2]
    end

    return options
end

# If variable expression is a symbol
# RV x ...
function extract_variable_id(expr::Symbol, options)
    if haskey(options, :id)
        return check_id_available(options[:id])
    else
        return guard_variable_id(:($(string(expr))))
    end
end

# If variable expression is an indexing expression
# RV x[i] ...
function extract_variable_id(expr::Expr, options)
    if haskey(options, :id)
        return check_id_available(options[:id])
    else
        argstr = map(arg -> :(string($arg)), @view expr.args[2:end])
        return guard_variable_id(:($(string(expr.args[1]) * "_") * $(reduce((current, item) -> :($current * "_" * $item), argstr))))
    end
end

# Fallback
function extract_variable_id(expr, options)
    return :(ForneyLab.generateId(Variable))
end


function check_id_available(expr)
    return :(!haskey(currentGraph().variables, $(expr)) ? $(expr) : error("Specified id is already assigned to another Variable"))
end

# Ensure that variable has a unique id in a current factor graph, generate a new one otherwise
function guard_variable_id(expr)
    idsymbol = :(Symbol($(expr)))
    return :(!haskey(currentGraph().variables, $idsymbol) ? $idsymbol : ForneyLab.generateId(Variable))
end