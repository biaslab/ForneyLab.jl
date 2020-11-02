export @ffg

using MacroTools: postwalk, rmlines, prettify, @capture

macro ffg(expr::Expr)
    return esc(postwalk(rmlines, generate_model(expr)))
end

function generate_model(expr::Expr)
    
    @capture(expr, (mname_(margs__) = body_) | (function mname_(margs__) body_ end)) || error("Model definition has to be a function.")

    body = postwalk(rewrite_expression, body)
    
    graph_sym = gensym(:factor_graph)

    result = quote
        function $mname($(margs...))
            $(graph_sym) = FactorGraph()
            $body
            return $(graph_sym)
        end
    end

    return result
end