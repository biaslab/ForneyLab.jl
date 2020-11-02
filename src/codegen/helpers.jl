# Extract options dictionary from expression
function get_options(expr::Expr)
    options = Dict{Symbol,Any}()
    
    expr = postwalk(expr) do x
        if @capture(x, lhs_ where {exoptions__})
            for option in exoptions
                @capture(option, key_ = value_)
                options[key] = value
            end
            return lhs
        end
        return x
    end

    if !isempty(options)
        return expr, options
    else
        return expr, nothing
    end
end

get_options(a::Any) = a, nothing