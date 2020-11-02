function rewrite_expression(expression::Expr)
    expr = if @capture(expression, var_ ~ rhs_)
        rewrite_tilde_expression(var, rhs)
    elseif @capture(expression, var_ = rhs_)
        rewrite_assign_expression(var, rhs)
    elseif is_for(expression)
        rewrite_for_block(expression)
    else
        expression
    end
    return expr
end

function rewrite_tilde_expression(var, rhs)
    
    (rhs, options) = get_options(rhs)
    @capture(rhs, pdist_(params__))
    
    var_id = extract_variable_id(var, options)
    
    # Build total expression
    return quote
        begin
            # Use existing Variable if it exists, otherwise create a new one
            $(var) = try
                $(var)
            catch _
                Variable(id = $(var_id))
            end

            # Create new variable if:
            #   - the existing object is not a Variable
            #   - the existing object is a Variable from another FactorGraph
            if (!isa($(var), Variable)
                || !haskey(currentGraph().variables, $(var).id)
                || currentGraph().variables[$(var).id] !== $(var))

                $(var) = Variable(id = $(var_id))
            end

            $(pdist)($(var), $(params...))
            $(var)
        end
    end
end

function rewrite_assign_expression(var, rhs)

    (rhs, options) = get_options(rhs)
    
    if options === nothing
        return quote $(var) = $(rhs) end
    end

    var_id = extract_variable_id(var, options)
    
    var_id_sym = gensym()

    return quote
        begin
            $(var) = $(rhs)
            $(var_id_sym) = $(var_id)
            if $(var_id_sym) != :auto
                # update id of newly created Variable
                currentGraph().variables[$(var_id_sym)] = $(var)
                delete!(currentGraph().variables, $(var).id)
                $(var).id = $(var_id_sym)
            end
            $(var)
        end
    end
end

# for loop
is_for(expr::Expr) = expr.head === :for
is_for(expr)       = false

function rewrite_for_block(def)
    body_block = def.args[2]
    
    for (i, expr) in enumerate(body_block.args)
        body_block.args[i] = rewrite_expression(expr)
    end
    
    return quote
        $(def)    
    end
end