function get_options(options_expr::Expr)
    options = Dict{Symbol,Any}()
    
    options_expr.head == :vect || return :(error("Incorrect options specification: options argument must be a vector expression"))
    
    for arg in options_expr.args
        arg isa Expr && arg.head == :(=) || return :(error("Incorrect options specification: options item must be an assignment expression"))
        options[arg.args[1]] = arg.args[2]
    end

    return options
end
