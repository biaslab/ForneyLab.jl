export CompositeNode, @composite

abstract type CompositeNode <: FactorNode end

"""
The `@composite` macro allows for defining custom (composite) nodes. Composite
nodes allow for implementating of custom update rules that may be computationally
more efficient or convenient. A composite node can be defined with or without an 
internal model. For detailed usage instructions we refer to the `composite_nodes` demo.
"""
macro composite(name::Symbol, exposed_vars::Expr, model::Expr)
    (exposed_vars.head == :tuple) || error("Exposed variables should be passed as Tuple")

    exposed_var_arguments = join(["$varname::Variable" for varname in exposed_vars.args], ", ")
    n_vars = length(exposed_vars.args)

    # Code for constructing interfaces
    interface_definitions = ""
    for idx = 1:n_vars
        varname = exposed_vars.args[idx]
        interface_definitions *= "self.i[:$varname] = self.interfaces[$idx] = ForneyLab.associate!(Interface(self), $varname)\n"
    end

    # Code for constructing exposed variables in inner graph
    exposed_var_definitions = ""
    terminal_definitions = ""
    for idx = 1:n_vars
        varname = exposed_vars.args[idx]
        exposed_var_definitions *= "local $varname = Variable(id=:$varname)\n"
        terminal_definitions *= "push!(self.terminals, Terminal($varname, self.interfaces[$idx], id=:$varname))\n"
        terminal_definitions *= "self.interface2terminal[self.interfaces[$idx]] = self.terminals[$idx]\n"
    end

    expr = parse("""
    mutable struct $name <: CompositeNode
        id::Symbol
        interfaces::Vector{Interface}
        i::Dict{Symbol, Interface}
        inner_graph::FactorGraph
        terminals::Vector{Terminal}
        interface2terminal::Dict{Interface,Terminal}

        function $name($exposed_var_arguments; id=ForneyLab.generateId($name))
            self = new(id, Array{Interface}(undef, $n_vars), Dict{Symbol,Interface}())
            outer_graph = currentGraph()
            ForneyLab.addNode!(outer_graph, self)
            $interface_definitions

            # Build internal graph
            self.inner_graph = FactorGraph()
            self.terminals = Terminal[]
            self.interface2terminal = Dict{Interface,Terminal}()
            let
                $exposed_var_definitions
                $model
                $terminal_definitions
            end
            setCurrentGraph(outer_graph)

            return self
        end
    end
    """)

    return esc(expr)
end