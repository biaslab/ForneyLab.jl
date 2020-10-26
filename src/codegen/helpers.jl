# Extract options dictionary from expression
function get_options(options_expr::Expr)
    options = Dict{Symbol,Any}()
    
    options_expr.head == :vect || return :(error("Incorrect options specification: options argument must be a vector expression"))
    
    for arg in options_expr.args
        arg isa Expr && arg.head == :(=) || return :(error("Incorrect options specification: options item must be an assignment expression"))
        options[arg.args[1]] = arg.args[2]
    end

    return options
end

# Check if options are defined for the expression
function tilde_options_are_defined(tilde_expr::Expr)
    try
        return tilde_expr.args[3].args[1] == :(∥)
    catch _
        return false
    end
end

function arrow_options_are_defined(arrow_expr::Expr)
    try
        return arrow_expr.args[3].args[1] == :(∥)
    catch _
        return false
    end
end

function assign_options_are_defined(assign_expr::Expr)
    try
        return assign_expr.args[2].args[1] == :(∥)
    catch _
        return false
    end
end
