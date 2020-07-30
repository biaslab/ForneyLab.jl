export Variable, @RV

"""
A `Variable` encompasses one or more edges in a `FactorGraph`.
"""
mutable struct Variable <: AbstractVariable
    id::Symbol
    edges::Vector{Edge}

    function Variable(;id=generateId(Variable))
        self = new(id, Vector{Edge}())
        addVariable!(currentGraph(), self)
        return self
    end
end

"""
`associate!(interface, variable)` associates `interface` with
`variable` by connecting `interface` to an `Edge` belonging to `variable`.
"""
function associate!(iface::Interface, var::Variable)
    if isempty(var.edges)
        # Make a new Edge
        Edge(var, iface)
    elseif var.edges[end].b == nothing
        # Connect to the loose end of an existing Edge
        connect!(var.edges[end], iface)
    else
        # Insert an equality constraint node
        equ_idx = Int((length(var.edges) - 1) / 2) + 1
        equ = Equality(id=Symbol("equ_$(var.id)_$(equ_idx)"))
        disconnected_iface = var.edges[end].b
        disconnect!(var.edges[end], disconnected_iface)
        connect!(var.edges[end], equ.interfaces[1])
        Edge(var, disconnected_iface, equ.interfaces[2])
        Edge(var, iface, equ.interfaces[3])
    end

    return iface
end

"""
Collect all edges corresponding with variable(s)
"""
edges(variable::Variable) = Set{Edge}(variable.edges)
edges(variables::Set{Variable}) = union(Set((edges(v) for v=variables))...)

Base.isless(v1::Variable, v2::Variable) = isless("$(v1.id)", "$(v2.id)")

"""
`@RV` provides a convenient way to add `Variable`s and `FactorNode`s to the graph.

Examples:

```
# Automatically create new Variable x, try to assign x.id = :x if this id is available
@RV x ~ GaussianMeanVariance(constant(0.0), constant(1.0))

# Explicitly specify the id of the Variable
@RV [id=:my_y] y ~ GaussianMeanVariance(constant(0.0), constant(1.0))

# Automatically assign z.id = :z if this id is not yet taken
@RV z = x + y

# Manual assignment
@RV [id=:my_u] u = x + y

# Just create a variable
@RV x
@RV [id=:my_x] x
```
"""
macro RV(options_expr::Expr, definition)
    # Parse options
    options_expr.head == :vect || return :(error("Incorrect use of @RV macro: options argument must be a vector expression"))
    definition isa Expr || definition isa Symbol || return :(error("Incorrect use of @RV macro: definition expression must be a valid expression or symbol"))
    options = Dict{Symbol, Any}()
    for arg in options_expr.args
        arg isa Expr && arg.head == :(=) || return :(error("Incorrect use of @RV macro: options item must be an assignment expression"))
        options[arg.args[1]] = arg.args[2]
    end

    # Parse RV definition expression
    # It can take three forms:
    # FORM 1: @RV x ~ Probdist(...)
    # FORM 2: @RV x = a + b
    # FORM 3: @RV x
    expr = if rv_isa_form1(definition)
        rv_form1(definition, definition.args[2], definition.args[3], options)
    elseif rv_isa_form2(definition)
        rv_form2(definition, definition.args[1], definition.args[2], options)
    elseif rv_isa_form3(definition)
        rv_form3(definition, definition, nothing, options)
    else
        :(error("Unsupported usage of @RV."))
    end

    return esc(expr)
end

macro RV(expr)
    # Pass empty options object
    return esc(:(@RV [] $(expr)))
end

# Parse RV definition expression

# FORM 1: @RV x ~ Probdist(...)
rv_isa_form1(expr::Expr) = expr.head === :call && expr.args[1] === :(~)
rv_isa_form1(expr)       = false

function rv_form1(def, target, node, options)
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

# FORM 2: @RV x = a + b
rv_isa_form2(expr::Expr) = expr.head === :(=)
rv_isa_form2(expr)       = false

function rv_form2(def, target, node, options)
    var_id = extract_variable_id(target, options)

    # Form 2 always creates a new Variable
    # Build complete expression
    var_id_sym = gensym()
    return quote
        begin
            $(def)
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

# FORM 3: @RV x
rv_isa_form3(expr::Symbol) = true
rv_isa_form3(expr::Expr)   = expr.head === :ref
rv_isa_form3(expr)         = false

function rv_form3(def, target, node, options)
    var_id = extract_variable_id(target, options)
    return quote
        $(target) = Variable(id = $(var_id))
    end
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
