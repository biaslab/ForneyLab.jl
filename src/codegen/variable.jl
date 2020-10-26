function rewrite_expression(expression::Expr)
    
    expr = if is_tilde(expression)
        rewrite_tilde_expression(expression)
    elseif is_arrow_assign(expression)
        rewrite_arrow_assign_expression(expression)
    elseif is_assign(expression)
        rewrite_assign_expression(expression)
    elseif is_for(expression)
        rewrite_for_block(expression)
    else
        expression
    end

    println("====================")
    dump(expr)
    println("********************")

    return expr
end


# Tilde expression: x ~ Probdist(...)
is_tilde(expr::Expr) = expr.head === :call && expr.args[1] === :(~)
is_tilde(expr)       = false

function rewrite_tilde_expression(def)
    if options_are_defined(def)
        options = get_options(def.args[3].args[3])
        node = def.args[3].args[2]
    else
        options = Dict{Symbol,Any}()
        node = def.args[3]
    end

    target = def.args[2]

    var_id = extract_variable_id(target, options)
    
    if isa(node.args[2], Expr) && (node.args[2].head == :parameters)
        node.args = vcat(node.args[1:2], [target], node.args[3:end])
    else
        node.args = vcat([node.args[1]; target], node.args[2:end])
    end
    
    # Build total expression
    return quote
        begin
            # Use existing Variable if it exists, otherwise create a new one
            $(target) = try
                $(target)
            catch _
                Variable(id = $(var_id))
            end

            # Create new variable if:
            #   - the existing object is not a Variable
            #   - the existing object is a Variable from another FactorGraph
            if (!isa($(target), Variable)
                || !haskey(currentGraph().variables, $(target).id)
                || currentGraph().variables[$(target).id] !== $(target))

                $(target) = Variable(id = $(var_id))
            end

            $(node)
            $(target)
        end
    end
end

# Arrow-style assignment: x ← a + b
is_arrow_assign(expr::Expr) = expr.head === :call && expr.args[1] === :(←)
is_arrow_assign(expr)       = false

function rewrite_arrow_assign_expression(def)
    if options_are_defined(def)
        options = get_options(def.args[3].args[3])
        rhs = def.args[3].args[2]
    else
        options = Dict{Symbol,Any}()
        rhs = def.args[3]
    end

    target = def.args[2]

    var_id = extract_variable_id(target, options)
    
    var_id_sym = gensym()

    return quote
        begin
            $(target) = $(rhs)
            $(var_id_sym) = $(var_id)
            if $(var_id_sym) != :auto
                # update id of newly created Variable
                currentGraph().variables[$(var_id_sym)] = $(target)
                delete!(currentGraph().variables, $(target).id)
                $(target).id = $(var_id_sym)
            end
            $(target)
        end
    end
end

# Regular assignment, possibly with options: x = a + b ∥ [id=:x]
is_assign(expr::Expr) = expr.head === :call && expr.args[1] === :(=)
is_assign(expr)       = false

function rewrite_assign_expression(def)
    return def
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