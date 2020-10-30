export @ffg

using MacroTools: postwalk, rmlines, prettify, @capture

macro ffg(ex::Expr)
    return esc(postwalk(rmlines, generate_model(ex)))
end

function generate_model(model_expr::Expr)
    program = postwalk(rmlines, model_expr)
    @assert program.head == :function
    
    model_signature = program.args[1]
    model_name, argument_names = analyze_signature(model_signature)
    
    model_definition = program.args[2]
    model_expr = build_model(model_definition)
    
    graph_sym = gensym(:factor_graph)

    result = quote
        function $model_name($(argument_names...))
            $(graph_sym) = FactorGraph()
            $model_expr
            return $(graph_sym)
        end
    end

    result = postwalk(rmlines, result)
    
    return result
end

function build_model(model_definition::Expr)
    for (i, expr) in enumerate(model_definition.args)
        model_definition.args[i] = rewrite_expression(expr)
    end
    
    return model_definition
end

function analyze_signature(args_expr)
    @assert args_expr.head == :call
    model_name = args_expr.args[1]
    
    if length(args_expr.args) > 1
        argument_names = args_expr.args[2:end]
    else
        argument_names = []
    end
    
    return model_name, argument_names
end


